## Pilot: how often do the two estimators fail, and how expensive is a
## replication?  Used only to choose the design grid and replication count for
## the main study.  Run from this directory:
##
##   cd dev/issue96 && Rscript 01-pilot.R

source('00-common.R')

pars <- make_pars()
Design <- createDesign(sample_size = c(50, 100, 200),
                       nitems = c(10, 20))

res <- runSimulation(Design, replications = 100, parallel = TRUE, ncores = 4,
                     generate = Generate,
                     analyse = list(FIML = Analyse.FIML, DWLS = Analyse.DWLS),
                     summarise = Summarise,
                     fixed_objects = pars,
                     control = list(allow_na = TRUE, allow_nan = TRUE))

print(as.data.frame(res)[, c('sample_size', 'nitems',
                             'FIML.converged', 'FIML.admissible', 'FIML.boundary',
                             'DWLS.converged', 'DWLS.postcheck', 'DWLS.admissible',
                             'DWLS.boundary', 'SIM_TIME')])
saveRDS(res, 'pilot.rds')
