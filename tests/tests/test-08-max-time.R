context('max_time')

test_that('max_time interruptions return partial results instead of crashing', {

    library(SimDesign)

    Design <- createDesign(N = 10)
    Generate <- function(condition, fixed_objects) rnorm(condition$N)
    AnaSlow <- function(condition, dat, fixed_objects){
        Sys.sleep(.2)
        c(m = mean(dat))
    }
    Summarise <- function(condition, results, fixed_objects)
        c(mu = mean(results[, "m"]))

    # a replication interrupted mid-run must be dropped like the pre-check
    # timed-out case (previously it produced a NULL that terminated with
    # "subscript out of bounds"); some replications complete before the
    # limit, so partial results are returned with a warning
    out <- suppressWarnings(
        runSimulation(Design, replications = 20, generate = Generate,
                      analyse = AnaSlow, summarise = Summarise, seed = 1,
                      parallel = FALSE, save = FALSE, verbose = FALSE,
                      control = list(max_time = 1)))
    expect_true(is.finite(out$mu))
    expect_lt(out$REPLICATIONS, 20)

    # when no replication completes at all (limit already expired before the
    # first replication starts) the condition must return a placeholder row
    # rather than crash the aggregation
    out2 <- suppressWarnings(
        runSimulation(Design, replications = 5, generate = Generate,
                      analyse = AnaSlow, summarise = Summarise, seed = 1,
                      parallel = FALSE, save = FALSE, verbose = FALSE,
                      control = list(max_time = .01)))
    expect_true("FATAL_TERMINATION" %in% colnames(out2))
    expect_true(grepl("max_time", out2$FATAL_TERMINATION[1]))

    # an interrupted replication must not overshoot the deadline by running
    # its retry unbounded (the caught interrupt clears the transient time
    # limit); with a 1s budget the 5s replication must be discarded rather
    # than completed on a second attempt. The busy loop supplies interpreter
    # checkpoints, which the time limit needs in order to interrupt (a single
    # long Sys.sleep() cannot be interrupted mid-sleep)
    AnaLong <- function(condition, dat, fixed_objects){
        for(i in 1:100){
            Sys.sleep(.05)
            invisible(sum(rnorm(10)))
        }
        c(m = mean(dat))
    }
    t0 <- proc.time()[3L]
    out3 <- suppressWarnings(
        runSimulation(Design, replications = 3, generate = Generate,
                      analyse = AnaLong, summarise = Summarise, seed = 1,
                      parallel = FALSE, save = FALSE, verbose = FALSE,
                      control = list(max_time = 1)))
    expect_lt(proc.time()[3L] - t0, 4)
    expect_true("FATAL_TERMINATION" %in% colnames(out3))

    # parallel workers must respect the wall-clock deadline; with the legacy
    # process-relative reference the freshly spawned workers computed inflated
    # budgets and evaluated every replication regardless of max_time
    out4 <- suppressWarnings(
        runSimulation(Design, replications = 30, generate = Generate,
                      analyse = AnaSlow, summarise = Summarise, seed = 1,
                      parallel = TRUE, ncores = 2, save = FALSE, verbose = FALSE,
                      control = list(max_time = 2)))
    expect_lt(out4$REPLICATIONS, 30)

})
