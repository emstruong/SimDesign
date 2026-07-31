## Validation run for philchalmers/SimDesign#96.
##
## 02/03 apply the redraw policy post hoc.  This script runs the redraw for
## real, using the vignette's own `stop()`-on-failure analysis functions, to
## confirm that (a) the post-hoc "listwise" numbers really are what redrawing
## produces, and (b) to record what the redrawing costs in discarded fits.
##
##   cd dev/issue96 && Rscript 04-redraw-run.R

source('00-common.R')

pars <- make_pars()

## Verbatim from vignettes/MultipleAnalyses.Rmd -----------------------------
Analyse.FIML.redraw <- function(condition, dat, fixed_objects) {
    mod <- mirt(dat, 1L, verbose = FALSE)
    if (!extract.mirt(mod, 'converged')) stop('mirt did not converge')
    cfs <- mirt::coef(mod, simplify = TRUE, digits = Inf)
    FIML_as <- cfs$items[, 1L] / 1.702
    c(as = unname(FIML_as))
}

Analyse.DWLS.redraw <- function(condition, dat, fixed_objects) {
    nitems <- condition$nitems
    lavmod <- paste0('F =~ ', paste0('NA*', colnames(dat)[1L], ' + '),
                     paste0(colnames(dat)[-1L], collapse = ' + '),
                     '\nF ~~ 1*F')
    lmod <- sem(lavmod, dat, ordered = colnames(dat))
    if (!lavInspect(lmod, 'converged')) stop('lavaan did not converge')
    cfs2 <- lavaan::coef(lmod)
    DWLS_alpha <- cfs2[1L:nitems]
    const <- sqrt(1 - DWLS_alpha^2)
    DWLS_as <- DWLS_alpha / const
    c(as = unname(DWLS_as))
}

## Same, but with the admissibility check requested in the issue -----------
Analyse.DWLS.redraw.adm <- function(condition, dat, fixed_objects) {
    nitems <- condition$nitems
    lavmod <- paste0('F =~ ', paste0('NA*', colnames(dat)[1L], ' + '),
                     paste0(colnames(dat)[-1L], collapse = ' + '),
                     '\nF ~~ 1*F')
    lmod <- sem(lavmod, dat, ordered = colnames(dat))
    if (!lavInspect(lmod, 'converged')) stop('lavaan did not converge')
    if (!lavaan::lavTech(lmod, 'post.check')) stop('lavaan solution not admissible')
    cfs2 <- lavaan::coef(lmod)
    DWLS_alpha <- cfs2[1L:nitems]
    const <- sqrt(1 - DWLS_alpha^2)
    DWLS_as <- DWLS_alpha / const
    c(as = unname(DWLS_as))
}

Summarise.recovery <- function(condition, results, fixed_objects) {
    nitems <- condition$nitems
    a <- fixed_objects[[ifelse(nitems == 10, 'ten', 'twenty')]]['a', ]
    pop <- c(a, a)
    c(bias = bias(results, pop), RMSE = RMSE(results, pop))
}

Design <- createDesign(sample_size = c(50, 100), nitems = c(10, 20))
REPLICATIONS <- 1000L

run_one <- function(dwls_fun, label) {
    t0 <- proc.time()[3L]
    res <- runSimulation(Design, replications = REPLICATIONS,
                         parallel = TRUE, ncores = 4,
                         generate = Generate,
                         analyse = list(FIML = Analyse.FIML.redraw,
                                        DWLS = dwls_fun),
                         summarise = Summarise.recovery,
                         fixed_objects = pars,
                         store_results = TRUE,
                         ## default max_errors = 50 halts a condition once the
                         ## redraw loop has burned through 50 tries; small-N
                         ## conditions legitimately need more than that
                         max_errors = 10000L,
                         seed = rep(20260801L, nrow(Design)) + seq_len(nrow(Design)),
                         verbose = TRUE)
    cat('\n---', label, '--- elapsed', round(proc.time()[3L] - t0, 1), 's\n')
    print(SimErrors(res))
    saveRDS(res, paste0('redraw-', label, '.rds'))
    res
}

res_vig <- run_one(Analyse.DWLS.redraw,     'vignette')
res_adm <- run_one(Analyse.DWLS.redraw.adm, 'admissible')
