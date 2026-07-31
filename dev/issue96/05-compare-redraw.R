## Checks that the post-hoc "listwise" policy in 03-analyse.R really is what
## an honest-to-goodness redraw produces, and prices the redraw in discarded
## model fits.
##
##   cd dev/issue96 && Rscript 05-compare-redraw.R

source('00-common.R')

pars <- make_pars()

to_loading <- function(a) a / sqrt(1 + a^2)

## Same reflection step as 03-analyse.R -- see the comment there.
reflect <- function(E) {
    s <- sign(rowSums(E, na.rm = TRUE))
    s[s == 0 | is.na(s)] <- 1
    E * s
}

summarise_est <- function(ests, pop) {
    bias_j <- colMeans(ests) - pop
    rmse_j <- sqrt(colMeans(sweep(ests, 2L, pop)^2))
    L <- to_loading(ests); popL <- to_loading(pop)
    c(n = nrow(ests), bias = mean(bias_j), RMSE = mean(rmse_j),
      bias_lam = mean(colMeans(L) - popL),
      RMSE_lam = mean(sqrt(colMeans(sweep(L, 2L, popL)^2))))
}

summarise_run <- function(file, label) {
    res <- readRDS(file)
    sr <- as.data.frame(SimResults(res))
    conds <- unique(sr[, c('sample_size', 'nitems')])
    out <- list()
    for (i in seq_len(nrow(conds))) {
        N <- conds$sample_size[i]; J <- conds$nitems[i]
        pop <- unname(pars[[ifelse(J == 10, 'ten', 'twenty')]]['a', ])
        d <- sr[sr$sample_size == N & sr$nitems == J, , drop = FALSE]
        for (est in c('FIML', 'DWLS')) {
            E <- reflect(as.matrix(d[, paste0(est, '.as', seq_len(J)),
                                     drop = FALSE]))
            s <- summarise_est(E, pop)
            out[[length(out) + 1L]] <- data.frame(
                run = label, sample_size = N, nitems = J, estimator = est,
                as.list(s), row.names = NULL)
        }
    }
    list(summary = do.call(rbind, out), errors = SimErrors(res),
         design = as.data.frame(res))
}

vig <- summarise_run('redraw-vignette.rds',   'redraw (vignette)')
adm <- summarise_run('redraw-admissible.rds', 'redraw (+post.check)')

## The post-hoc listwise numbers from the keep-everything run
pc <- readRDS('policy-comparison.rds')
lw <- pc$tab[pc$tab$policy == 'listwise',
             c('sample_size', 'nitems', 'estimator', 'n', 'bias', 'RMSE',
               'bias_lam', 'RMSE_lam')]
lw <- cbind(run = 'post-hoc listwise', lw)

keep <- pc$tab[pc$tab$policy == 'marginal',
               c('sample_size', 'nitems', 'estimator', 'n', 'bias', 'RMSE',
                 'bias_lam', 'RMSE_lam')]
keep <- cbind(run = 'post-hoc marginal (no redraw)', keep)

cmp <- rbind(vig$summary, adm$summary, lw, keep)
cmp <- cmp[cmp$sample_size %in% c(50, 100), ]
cmp <- cmp[order(cmp$nitems, cmp$sample_size, cmp$estimator, cmp$run), ]

options(digits = 4)
cat('\n==== actual redraw runs vs the post-hoc policies ====\n')
print(cmp, row.names = FALSE)

cat('\n==== errors consumed by the redraw loop (vignette policy) ====\n')
print(vig$errors)
cat('\n==== errors consumed by the redraw loop (+post.check) ====\n')
print(adm$errors)

cat('\n==== wall time ====\n')
print(data.frame(run = 'redraw (vignette)',
                 vig$design[, c('sample_size', 'nitems', 'SIM_TIME')]),
      row.names = FALSE)
print(data.frame(run = 'redraw (+post.check)',
                 adm$design[, c('sample_size', 'nitems', 'SIM_TIME')]),
      row.names = FALSE)

write.csv(cmp, 'redraw-validation.csv', row.names = FALSE)
