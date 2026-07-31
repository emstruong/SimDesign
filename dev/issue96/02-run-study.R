## Main study for philchalmers/SimDesign#96.
##
## One expensive "keep everything" run.  Because every replication is retained
## along with per-estimator convergence/admissibility flags, any redraw policy
## can be applied *post hoc* to exactly the same Monte Carlo draws.  That makes
## the policy comparison paired: the difference between policies cannot be
## attributed to Monte Carlo noise between two independent runs.
##
##   cd dev/issue96 && Rscript 02-run-study.R

source('00-common.R')

pars <- make_pars()

## The vignette itself uses sample_size = 250, 500, 1000.  The grid is extended
## downwards because that is where the estimators actually fail.
Design <- createDesign(sample_size = c(50, 100, 250, 500, 1000),
                       nitems = c(10, 20))

REPLICATIONS <- 2000L

res <- runSimulation(Design, replications = REPLICATIONS,
                     parallel = TRUE, ncores = 4,
                     generate = Generate,
                     analyse = list(FIML = Analyse.FIML, DWLS = Analyse.DWLS),
                     summarise = Summarise,
                     fixed_objects = pars,
                     store_results = TRUE,
                     seed = rep(20260731L, nrow(Design)) + seq_len(nrow(Design)),
                     ## allow_na/allow_nan are what switch the redraw *off*:
                     ## an unusable estimate is returned as NA and kept.
                     control = list(allow_na = TRUE, allow_nan = TRUE),
                     verbose = TRUE)

saveRDS(res, 'keep-run.rds')
saveRDS(SimResults(res), 'keep-results.rds')

print(as.data.frame(res)[, c('sample_size', 'nitems',
                             'FIML.converged', 'FIML.admissible', 'FIML.boundary',
                             'DWLS.converged', 'DWLS.postcheck', 'DWLS.admissible',
                             'DWLS.boundary', 'SIM_TIME')])
