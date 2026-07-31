context('error-tally')

test_that('fuzzy_reduce does not double count', {

    # Three distinct error messages, where the short DWLS-style message also
    # occurs inside the longer combined message.  Grouping the short one must
    # not additionally absorb the combined column, which was already reported
    # as its own group.
    df <- dplyr::tibble(a = 118L, b = 7L, c = 250L)
    colnames(df) <- c(
        paste0("ERROR:  2 INDEPENDENT ERRORS THROWN:   FIML.ERROR:   ",
               "mirt did not converge\n  DWLS.ERROR:   ",
               "lavaan solution not admissible\n\n"),
        "ERROR:  DWLS.ERROR:   lavaan solution not admissible\n\n",
        "ERROR:  FIML.ERROR:   mirt did not converge\n\n")

    out <- SimDesign:::fuzzy_reduce(df)
    expect_equal(sum(unlist(out)), sum(unlist(df)))
    expect_equal(unname(unlist(out)), c(118L, 7L, 250L))

    # a single column is returned untouched
    one <- dplyr::tibble(x = 5L)
    colnames(one) <- "ERROR:  boom\n\n"
    expect_equal(sum(unlist(SimDesign:::fuzzy_reduce(one))), 5L)

    # messages far apart are kept separate and summed correctly
    far <- dplyr::tibble(a = 3L, b = 4L)
    colnames(far) <- c("ERROR:  completely unrelated failure\n\n",
                       "ERROR:  something else entirely different here\n\n")
    expect_equal(sum(unlist(SimDesign:::fuzzy_reduce(far))), 7L)
})

test_that('SimErrors tally matches the number of re-draws', {

    Design <- createDesign(dummy = 1)
    Generate <- function(condition, fixed_objects) rnorm(1)
    # A fails for x > 1; B fails for x > 1.5 (both fail) or x < -0.5 (B alone)
    Analyse.A <- function(condition, dat, fixed_objects) {
        if(dat > 1) stop('A failed')
        c(a = dat)
    }
    Analyse.B <- function(condition, dat, fixed_objects) {
        if(dat > 1.5 || dat < -0.5) stop('B failed')
        c(b = dat)
    }
    Summarise <- function(condition, results, fixed_objects) c(m = mean(results[,1]))

    res <- runSimulation(Design, replications = 200, generate = Generate,
                         analyse = list(A = Analyse.A, B = Analyse.B),
                         summarise = Summarise, max_errors = 100000L,
                         seed = 4321L, verbose = FALSE)

    # the ERRORS column counts re-draws; the per-message table must agree
    tally <- sum(as.numeric(SimErrors(res)[1L, -1L, drop = FALSE]))
    expect_equal(tally, res$ERRORS)

    # fuzzy = FALSE reports every message separately and must also agree
    exact <- sum(as.numeric(SimErrors(res, fuzzy = FALSE)[1L, -1L, drop = FALSE]))
    expect_equal(exact, res$ERRORS)
})

test_that('NA/NaN re-draw messages list every offending element', {

    Design <- createDesign(dummy = 1)
    Generate <- function(condition, fixed_objects) rnorm(1)
    Analyse <- function(condition, dat, fixed_objects) c(a = NaN, b = NaN, c = 1)
    Summarise <- function(condition, results, fixed_objects) c(m = 1)

    res <- try(runSimulation(Design, replications = 2, generate = Generate,
                             analyse = Analyse, summarise = Summarise,
                             max_errors = 3L, verbose = FALSE), silent = TRUE)
    msg <- as.character(res)[1L]
    expect_true(grepl('redrawing: a, b', msg, fixed = TRUE))
})
