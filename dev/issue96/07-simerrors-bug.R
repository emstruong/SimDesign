## Self-contained reproducer for the two SimErrors() message-table problems
## found while validating the re-draw runs for issue #96.
##
##   Rscript 07-simerrors-bug.R
##
## Both stem from fuzzy_reduce() grouping messages with agrepl(), whose default
## max.distance is 10% of the pattern length and which matches a pattern
## occurring *inside* a longer string.

suppressMessages(library(SimDesign))

## The old implementation, for side-by-side comparison ----------------------
old_fuzzy_reduce <- function(df){
    nms <- colnames(df); matched <- logical(length(nms)); udf <- df[, 0]
    for(i in seq_along(nms)){
        if(matched[i]) next
        udf <- cbind(udf, df[, i])
        tm <- agrepl(nms[i], nms)
        udf[, ncol(udf)] <- rowSums(df[, tm, drop = FALSE], na.rm = TRUE)
        matched <- matched | tm
    }
    udf
}

## Problem 1: over-counting ------------------------------------------------
## The short DWLS message occurs inside the combined message, so the combined
## column is reported as its own group and summed into the DWLS group again.

df <- dplyr::tibble(both = 118L, dwls_only = 7L, mirt_only = 250L)
colnames(df) <- c(
    paste0("ERROR:  2 INDEPENDENT ERRORS THROWN:   FIML.ERROR:   ",
           "mirt did not converge\n  DWLS.ERROR:   ",
           "lavaan solution not admissible\n\n"),
    "ERROR:  DWLS.ERROR:   lavaan solution not admissible\n\n",
    "ERROR:  FIML.ERROR:   mirt did not converge\n\n")

cat('=== Problem 1: totals exceed the number of re-draws ===\n')
cat('true counts        :', unlist(df), ' total', sum(unlist(df)), '\n')
cat('old fuzzy_reduce   :', unlist(old_fuzzy_reduce(df)),
    ' total', sum(unlist(old_fuzzy_reduce(df))), '\n')
cat('fixed fuzzy_reduce :', unlist(SimDesign:::fuzzy_reduce(df)),
    ' total', sum(unlist(SimDesign:::fuzzy_reduce(df))), '\n')

## Problem 2: mis-attribution ----------------------------------------------
## Two unrelated failures whose messages differ by a handful of characters are
## inside agrep's default tolerance, so they are merged under the first name.

cat('\n=== Problem 2: distinct errors merged under the first name ===\n')
a <- "ERROR:  A.ERROR:   AAA failed\n\n"
b <- "ERROR:  B.ERROR:   BBB failed\n\n"
cat('agrepl("A.ERROR: AAA failed", "B.ERROR: BBB failed") =', agrepl(a, b), '\n')

Design <- createDesign(dummy = 1)
Generate <- function(condition, fixed_objects) rnorm(1)
## disjoint failure regions: A and B can never fail on the same attempt
Analyse.A <- function(condition, dat, fixed_objects) {
    if(dat >  1.0) stop('AAA failed'); c(a = dat)
}
Analyse.B <- function(condition, dat, fixed_objects) {
    if(dat < -0.5) stop('BBB failed'); c(b = dat)
}
Summarise <- function(condition, results, fixed_objects) c(m = mean(results[, 1]))

res <- runSimulation(Design, replications = 500, generate = Generate,
                     analyse = list(A = Analyse.A, B = Analyse.B),
                     summarise = Summarise, max_errors = 100000L,
                     seed = 1L, verbose = FALSE)

cat('\nSimErrors() default (fuzzy = TRUE) -- one column, wrong attribution:\n')
print(t(SimErrors(res)))
cat('\nSimErrors(fuzzy = FALSE) -- the actual per-function breakdown:\n')
print(t(SimErrors(res, fuzzy = FALSE)))
cat('\nunderlying ERROR_msg attribute (always correct):\n')
print(unlist(attr(res, 'ERROR_msg')[1, ]))
