context('load_balancing')

test_that('load balanced dispatches are reproducible', {

    library(SimDesign)

    Design <- createDesign(N = c(10, 20))

    Generate <- function(condition, fixed_objects) rnorm(condition$N)
    Analyse <- function(condition, dat, fixed_objects) c(mean = mean(dat))
    Summarise <- function(condition, results, fixed_objects)
        c(mu = mean(results[, "mean"]))

    # same seed must give identical results regardless of the number of
    # cores or the backend used, as RNG streams are tied to replications
    # rather than to worker nodes
    res_mirai2 <- runSimulation(Design, replications = 6, generate = Generate,
                                analyse = Analyse, summarise = Summarise,
                                seed = c(123, 456), parallel = TRUE, ncores = 2,
                                save = FALSE, verbose = FALSE)
    res_mirai3 <- runSimulation(Design, replications = 6, generate = Generate,
                                analyse = Analyse, summarise = Summarise,
                                seed = c(123, 456), parallel = TRUE, ncores = 3,
                                save = FALSE, verbose = FALSE)
    expect_equal(res_mirai2$mu, res_mirai3$mu)
    res_psock <- runSimulation(Design, replications = 6, generate = Generate,
                               analyse = Analyse, summarise = Summarise,
                               seed = c(123, 456), parallel = TRUE, ncores = 2,
                               save = FALSE, verbose = FALSE,
                               control = list(use_mirai = FALSE))
    expect_equal(res_mirai2$mu, res_psock$mu)

    # legacy static scheduling remains available and reproducible
    res_static <- runSimulation(Design, replications = 6, generate = Generate,
                                analyse = Analyse, summarise = Summarise,
                                seed = c(123, 456), parallel = TRUE, ncores = 2,
                                save = FALSE, verbose = FALSE,
                                control = list(use_load_balancing = FALSE))
    res_static2 <- runSimulation(Design, replications = 6, generate = Generate,
                                 analyse = Analyse, summarise = Summarise,
                                 seed = c(123, 456), parallel = TRUE, ncores = 2,
                                 save = FALSE, verbose = FALSE,
                                 control = list(use_load_balancing = FALSE))
    expect_equal(res_static$mu, res_static2$mu)

    # fatal terminations are still collected under dynamic scheduling
    AnalyseBad <- function(condition, dat, fixed_objects) stop('always fails')
    expect_warning(res_bad <- runSimulation(Design[1L, , drop=FALSE], replications = 6,
                             generate = Generate, analyse = AnalyseBad,
                             summarise = Summarise, seed = 42, max_errors = 2,
                             parallel = TRUE, ncores = 2,
                             save = FALSE, verbose = FALSE,
                             control = list(stop_on_fatal = FALSE)),
                   "fatally terminated")
    expect_true("FATAL_TERMINATION" %in% colnames(res_bad))

    # runArraySimulation() sub-stream seeds are also invariant to ncores.
    # Note that functions are placed in the global environment here, as
    # the array + parallel combination discovers exportable functions from
    # the caller but resolves them lexically (both scheduling approaches)
    tmpdir <- tempfile()
    dir.create(tmpdir)
    owd <- setwd(tmpdir)
    on.exit(setwd(owd), add = TRUE)
    assign('Generate.arr', Generate, envir = globalenv())
    assign('Analyse.arr', Analyse, envir = globalenv())
    assign('Summarise.arr', Summarise, envir = globalenv())
    on.exit(rm('Generate.arr', 'Analyse.arr', 'Summarise.arr',
               envir = globalenv()), add = TRUE)
    arr2 <- with(globalenv(),
        runArraySimulation(createDesign(N = c(10, 20)), replications = 6,
                           generate = Generate.arr,
                           analyse = Analyse.arr, summarise = Summarise.arr,
                           iseed = 554184288, arrayID = 2, filename = 'arr2',
                           parallel = TRUE, ncores = 2, verbose = FALSE))
    arr3 <- with(globalenv(),
        runArraySimulation(createDesign(N = c(10, 20)), replications = 6,
                           generate = Generate.arr,
                           analyse = Analyse.arr, summarise = Summarise.arr,
                           iseed = 554184288, arrayID = 2, filename = 'arr3',
                           parallel = TRUE, ncores = 3, verbose = FALSE))
    expect_equal(arr2$mu, arr3$mu)
    SimClean(dir())

})
