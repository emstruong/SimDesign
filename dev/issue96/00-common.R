## Common Generate/Analyse/Summarise definitions for the issue #96 study.
##
## The simulation is the parameter-recovery example from
## vignettes/MultipleAnalyses.Rmd: a unidimensional normal-ogive IRT model is
## fitted with marginal ML (`mirt`, "FIML") and with an ordinal CFA using
## diagonally weighted least squares (`lavaan`, "DWLS").  Both estimators
## recover the same target: the normal-ogive item slopes `a`.
##
## The study adds small-N conditions, and records *per replication* whether
## each estimator converged and whether its solution was admissible.  Those
## flags let the redraw policy be applied post hoc, so a single expensive run
## supports every policy comparison.

library(SimDesign)
library(mirt)
library(lavaan)

## ---- population values --------------------------------------------------
## Identical to the vignette (set.seed(1) draws).
make_pars <- function() {
    set.seed(1)
    pars_10 <- rbind(a = round(rlnorm(10, .3, .5) / 1.702, 2),
                     d = round(rnorm(10, 0, .5) / 1.702, 2))
    pars_20 <- rbind(a = round(rlnorm(20, .3, .5) / 1.702, 2),
                     d = round(rnorm(20, 0, .5) / 1.702, 2))
    list(ten = pars_10, twenty = pars_20)
}

P_ogive <- function(a, d, Theta) pnorm(a * Theta + d)

## `a` estimates larger than this are treated as boundary/degenerate solutions
## in the sensitivity analysis.  In the normal-ogive metric a = 5 corresponds
## to a standardised loading of a/sqrt(1 + a^2) = .98.
BOUNDARY_A <- 5

## ---- Generate -----------------------------------------------------------
Generate <- function(condition, fixed_objects) {
    N <- condition$sample_size
    nitems <- condition$nitems
    nitems_name <- ifelse(nitems == 10, 'ten', 'twenty')

    a <- fixed_objects[[nitems_name]]['a', ]
    d <- fixed_objects[[nitems_name]]['d', ]

    Theta <- rnorm(N)
    dat <- matrix(0L, N, nitems)
    for (j in 1:nitems) {
        p <- P_ogive(a[j], d[j], Theta)
        dat[, j] <- rbinom(N, 1L, p)
    }
    colnames(dat) <- paste0('item_', 1:nitems)
    as.data.frame(dat)
}

## ---- Analyse ------------------------------------------------------------
## Both analysis functions return
##   as1 ... as<nitems>  the recovered normal-ogive slopes (NA if unusable)
##   converged           1/0, the estimator's own convergence flag
##   admissible          1/0, converged *and* the solution passes the
##                       estimator's admissibility check
##   boundary            1/0, any |a-hat| >= BOUNDARY_A (sensitivity analysis)
## Nothing is ever `stop()`ped, so SimDesign never redraws: every replication
## is retained together with the flags describing what happened to it.

Analyse.FIML <- function(condition, dat, fixed_objects) {
    nitems <- condition$nitems
    as <- rep(NA_real_, nitems)
    converged <- 0
    boundary <- 0

    mod <- try(mirt(dat, 1L, verbose = FALSE), silent = TRUE)
    if (!is(mod, 'try-error')) {
        converged <- as.numeric(extract.mirt(mod, 'converged'))
        cfs <- try(mirt::coef(mod, simplify = TRUE, digits = Inf), silent = TRUE)
        if (!is(cfs, 'try-error')) {
            FIML_as <- unname(cfs$items[, 1L] / 1.702)
            if (length(FIML_as) == nitems) as <- FIML_as
        }
    }
    ## The IRT parameterisation has no Heywood case: the implied standardised
    ## loading a/sqrt(1 + a^2) is bounded by 1 for any finite a.  The
    ## degenerate solution is instead a -> Inf (quasi-separation), so
    ## "admissible" for FIML is convergence plus finite, non-boundary slopes.
    finite_ok <- all(is.finite(as))
    boundary <- as.numeric(finite_ok && any(abs(as) >= BOUNDARY_A))
    admissible <- as.numeric(converged == 1 && finite_ok)

    c(as = as, converged = converged, admissible = admissible, boundary = boundary)
}

Analyse.DWLS <- function(condition, dat, fixed_objects) {
    nitems <- condition$nitems
    as <- rep(NA_real_, nitems)
    converged <- 0
    postcheck <- 0
    boundary <- 0

    lavmod <- paste0('F =~ ', paste0('NA*', colnames(dat)[1L], ' + '),
                     paste0(colnames(dat)[-1L], collapse = ' + '),
                     '\nF ~~ 1*F')
    lmod <- try(suppressWarnings(sem(lavmod, dat, ordered = colnames(dat))),
                silent = TRUE)

    if (!is(lmod, 'try-error')) {
        converged <- as.numeric(isTRUE(lavInspect(lmod, 'converged')))
        ## `post.check` is lavaan's own admissibility test (no negative
        ## variances, no |correlation| > 1, positive definite implied
        ## matrices).  It is *independent* of the convergence flag: lavaan
        ## routinely converges on an inadmissible solution.
        pc <- try(lavaan::lavTech(lmod, 'post.check'), silent = TRUE)
        postcheck <- if (is(pc, 'try-error')) 0 else as.numeric(isTRUE(pc))

        cfs2 <- try(lavaan::coef(lmod), silent = TRUE)
        if (!is(cfs2, 'try-error') && length(cfs2) >= nitems) {
            DWLS_alpha <- cfs2[1L:nitems]
            ## sqrt(1 - alpha^2) is NaN for a Heywood loading; keep the NaN so
            ## the failure is visible rather than silently dropped.
            const <- suppressWarnings(sqrt(1 - DWLS_alpha^2))
            as <- unname(DWLS_alpha / const)
        }
    }
    finite_ok <- all(is.finite(as))
    boundary <- as.numeric(finite_ok && any(abs(as) >= BOUNDARY_A))
    ## An admissible DWLS solution must converge, pass post.check, and yield
    ## a usable slope transformation (|loading| < 1).
    admissible <- as.numeric(converged == 1 && postcheck == 1 && finite_ok)

    c(as = as, converged = converged, admissible = admissible,
      boundary = boundary, postcheck = postcheck)
}

## ---- Summarise ----------------------------------------------------------
## Deliberately minimal: all policy comparisons are done post hoc from the
## stored per-replication results.  Only the flag rates are summarised here so
## that the printed SimDesign object is informative on its own.
Summarise <- function(condition, results, fixed_objects) {
    nms <- colnames(results)
    flag <- nms[grepl('(converged|admissible|boundary|postcheck)$', nms)]
    colMeans(results[, flag, drop = FALSE], na.rm = TRUE)
}
