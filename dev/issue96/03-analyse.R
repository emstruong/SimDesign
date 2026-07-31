## Post-hoc policy comparison for philchalmers/SimDesign#96.
##
## Reads the per-replication results from 02-run-study.R and evaluates three
## policies on *the same* Monte Carlo draws:
##
##   none      every replication the estimator returns a finite value for,
##             regardless of convergence or admissibility
##   marginal  drop only the replications where *this* estimator failed
##             (what you get when the simulation does NOT redraw)
##   listwise  drop the replication whenever *either* estimator failed
##             (what redrawing produces: a redraw discards the other
##             estimator's perfectly good fit along with the bad one)
##
## Because all three are computed from one set of draws, the marginal-vs-
## listwise contrast is paired and cannot be explained by run-to-run noise.
##
##   cd dev/issue96 && Rscript 03-analyse.R

source('00-common.R')
set.seed(20260731)

pars <- make_pars()
sr <- as.data.frame(readRDS('keep-results.rds'))
B_BOOT <- 1000L

## Two scales for the same estimand.  The vignette reports the normal-ogive
## slope `a`, which is unbounded and explodes as the implied loading
## approaches 1, so a handful of near-boundary replications dominate its RMSE.
## The standardised loading lambda = a / sqrt(1 + a^2) is a monotone
## reparameterisation of the identical quantity but is bounded by 1, so it
## shows the policy effect without the heavy tail.  Both are reported.
to_loading <- function(a) a / sqrt(1 + a^2)

## `F ~~ 1*F` fixes the factor's variance but not its orientation, so the
## reflected solution fits identically and lavaan returns it a nontrivial
## share of the time once N is small.  That is an identification artifact, not
## a convergence or admissibility failure -- lavaan reports such solutions as
## converged and post.check-clean -- and it is unrelated to the redraw
## question, so it is resolved here before anything is compared to the (all
## positive) population slopes.  `mirt` orients the factor internally and
## never needs this.  Doing it here rather than inside Analyse.DWLS keeps
## 00-common.R's analysis functions faithful to the vignette's code.
reflect <- function(E) {
    s <- sign(rowSums(E, na.rm = TRUE))
    s[s == 0 | is.na(s)] <- 1
    E * s
}

est_cols <- function(dat, est, nitems)
    reflect(as.matrix(dat[, paste0(est, '.as', seq_len(nitems)), drop = FALSE]))

## Policy masks -------------------------------------------------------------
## `usable` is a hard requirement under every policy: a non-finite estimate
## cannot enter a mean or an RMSE regardless of the redraw decision.
policy_mask <- function(dat, est, nitems, policy) {
    other <- setdiff(c('FIML', 'DWLS'), est)
    E <- est_cols(dat, est, nitems)
    usable <- apply(is.finite(E), 1L, all)
    switch(policy,
           none     = usable,
           marginal = usable & dat[[paste0(est, '.admissible')]] == 1,
           listwise = usable &
               dat[[paste0(est,   '.admissible')]] == 1 &
               dat[[paste0(other, '.admissible')]] == 1,
           stop('unknown policy'))
}

summarise_est <- function(ests, pop) {
    if (nrow(ests) == 0L)
        return(c(n = 0, bias = NA, abs_bias = NA, RMSE = NA,
                 bias_lam = NA, abs_bias_lam = NA, RMSE_lam = NA))
    bias_j <- colMeans(ests) - pop
    rmse_j <- sqrt(colMeans(sweep(ests, 2L, pop)^2))
    L <- to_loading(ests)
    popL <- to_loading(pop)
    biasL_j <- colMeans(L) - popL
    rmseL_j <- sqrt(colMeans(sweep(L, 2L, popL)^2))
    c(n = nrow(ests), bias = mean(bias_j), abs_bias = mean(abs(bias_j)),
      RMSE = mean(rmse_j), bias_lam = mean(biasL_j),
      abs_bias_lam = mean(abs(biasL_j)), RMSE_lam = mean(rmseL_j))
}

evaluate <- function(dat, est, nitems, pop, policy) {
    m <- policy_mask(dat, est, nitems, policy)
    summarise_est(est_cols(dat, est, nitems)[m, , drop = FALSE], pop)
}

## Bootstrap over replications gives the standard error of the *paired*
## marginal-minus-listwise difference, i.e. the Monte Carlo error the
## difference has to clear to mean anything.
boot_diff <- function(dat, est, nitems, pop, B = B_BOOT) {
    keys <- c('bias', 'RMSE', 'bias_lam', 'RMSE_lam')
    n <- nrow(dat)
    out <- matrix(NA_real_, B, length(keys), dimnames = list(NULL, keys))
    for (b in seq_len(B)) {
        d <- dat[sample.int(n, n, replace = TRUE), , drop = FALSE]
        mg <- evaluate(d, est, nitems, pop, 'marginal')
        lw <- evaluate(d, est, nitems, pop, 'listwise')
        out[b, ] <- mg[keys] - lw[keys]
    }
    apply(out, 2L, stats::sd, na.rm = TRUE)
}

conditions <- unique(sr[, c('sample_size', 'nitems')])
rows <- diffs <- fails <- list()

for (i in seq_len(nrow(conditions))) {
    N <- conditions$sample_size[i]
    J <- conditions$nitems[i]
    pop <- unname(pars[[ifelse(J == 10, 'ten', 'twenty')]]['a', ])
    dat <- sr[sr$sample_size == N & sr$nitems == J, , drop = FALSE]

    ## Which estimator is responsible for the discarded replications?
    fok <- dat$FIML.admissible == 1
    dok <- dat$DWLS.admissible == 1
    rawD <- as.matrix(dat[, paste0('DWLS.as', seq_len(J)), drop = FALSE])
    fails[[length(fails) + 1L]] <- data.frame(
        sample_size = N, nitems = J, reps = nrow(dat),
        DWLS_reflected = mean(rowSums(rawD, na.rm = TRUE) < 0),
        both_ok = mean(fok & dok), FIML_only_fails = mean(!fok & dok),
        DWLS_only_fails = mean(fok & !dok), both_fail = mean(!fok & !dok),
        FIML_fail = mean(!fok), DWLS_fail = mean(!dok),
        ## the share of DWLS's usable fits that a redraw throws away
        DWLS_collateral = if (any(dok)) mean(!fok[dok]) else NA_real_,
        FIML_collateral = if (any(fok)) mean(!dok[fok]) else NA_real_,
        row.names = NULL)

    for (est in c('FIML', 'DWLS')) {
        for (policy in c('none', 'marginal', 'listwise')) {
            s <- evaluate(dat, est, J, pop, policy)
            rows[[length(rows) + 1L]] <- data.frame(
                sample_size = N, nitems = J, estimator = est, policy = policy,
                as.list(s), row.names = NULL)
        }
        se <- boot_diff(dat, est, J, pop)
        mg <- evaluate(dat, est, J, pop, 'marginal')
        lw <- evaluate(dat, est, J, pop, 'listwise')
        diffs[[length(diffs) + 1L]] <- data.frame(
            sample_size = N, nitems = J, estimator = est,
            n_marginal = mg[['n']], n_listwise = lw[['n']],
            d_bias = mg[['bias']] - lw[['bias']], se_bias = se[['bias']],
            d_RMSE = mg[['RMSE']] - lw[['RMSE']], se_RMSE = se[['RMSE']],
            d_bias_lam = mg[['bias_lam']] - lw[['bias_lam']],
            se_bias_lam = se[['bias_lam']],
            d_RMSE_lam = mg[['RMSE_lam']] - lw[['RMSE_lam']],
            se_RMSE_lam = se[['RMSE_lam']], row.names = NULL)
    }
}

tab <- do.call(rbind, rows)
dif <- do.call(rbind, diffs)
fail <- do.call(rbind, fails)
dif$z_bias_lam <- dif$d_bias_lam / dif$se_bias_lam
dif$z_RMSE_lam <- dif$d_RMSE_lam / dif$se_RMSE_lam
dif$z_bias <- dif$d_bias / dif$se_bias

## Does the *conclusion* change?  The vignette compares two estimators, so the
## conclusion of interest is the gap between them, not either one alone.
gap <- do.call(rbind, lapply(
    split(tab, list(tab$sample_size, tab$nitems, tab$policy), drop = TRUE),
    function(d) data.frame(
        sample_size = d$sample_size[1L], nitems = d$nitems[1L],
        policy = d$policy[1L],
        RMSE_lam_FIML = d$RMSE_lam[d$estimator == 'FIML'],
        RMSE_lam_DWLS = d$RMSE_lam[d$estimator == 'DWLS'],
        bias_lam_FIML = d$bias_lam[d$estimator == 'FIML'],
        bias_lam_DWLS = d$bias_lam[d$estimator == 'DWLS'],
        row.names = NULL)))
gap$RMSE_gap <- gap$RMSE_lam_FIML - gap$RMSE_lam_DWLS
gap$winner <- ifelse(gap$RMSE_gap < 0, 'FIML', 'DWLS')
gap <- gap[order(gap$nitems, gap$sample_size, gap$policy), ]

ord <- function(d) {
    key <- if (is.null(d$estimator)) rep('', nrow(d)) else d$estimator
    d[order(d$nitems, d$sample_size, key), ]
}

options(digits = 4)
cat('\n==== failure rates and who causes the discards ====\n')
print(ord(fail), row.names = FALSE)
cat('\n==== per-estimator summaries by policy ====\n')
print(tab[order(tab$nitems, tab$sample_size, tab$estimator, tab$policy), ],
      row.names = FALSE)
cat('\n==== marginal (no redraw) minus listwise (redraw), bootstrap SE ====\n')
print(ord(dif), row.names = FALSE)
cat('\n==== estimator comparison by policy (loading scale) ====\n')
print(gap, row.names = FALSE)

saveRDS(list(tab = tab, dif = dif, gap = gap, fail = fail),
        'policy-comparison.rds')
write.csv(tab,  'policy-summaries.csv',   row.names = FALSE)
write.csv(dif,  'policy-differences.csv', row.names = FALSE)
write.csv(gap,  'policy-rmse-gap.csv',    row.names = FALSE)
write.csv(fail, 'failure-rates.csv',      row.names = FALSE)
