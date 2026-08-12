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

})
