.SIMDENV <- new.env(parent=emptyenv())

# Return a character vector of functions defined in the calling stack frames,
# with the discovered function objects captured in an 'envir' attribute.
# Discovery walks the *call stack* (parent.frame(level:2)), yet the export
# step historically resolved the returned names through a single frame's
# *lexical* chain (envir = parent.frame(1L) at the clusterExport() call),
# which errors with "object not found" when runSimulation() is reached
# through a wrapper such as runArraySimulation() and the user's functions
# live outside .GlobalEnv (e.g., defined within a function body or a
# testthat block). Capturing the objects at discovery time keeps the
# lookup and the export consistent
parent_env_fun <- function(level=2){
    ret <- NULL
    objs <- new.env(parent = emptyenv())
    for(lev in level:2){
        nms <- ls(envir = parent.frame(lev))
        is_fun <- sapply(nms, function(x, envir) is.function(get(x, envir=envir)),
                         envir = parent.frame(lev))
        if(any(is_fun)){
            ret <- c(ret, nms[is_fun])
            # frames closer to the runSimulation() call take precedence
            # over more distant ones when names are duplicated, matching
            # the historical get() resolution order
            for(nm in nms[is_fun])
                assign(nm, get(nm, envir = parent.frame(lev)), envir = objs)
        }
    }
    if(!is.null(ret)){
        ret <- unique(ret)
        attr(ret, 'envir') <- objs
    }
    ret
}

unique_filename <- function(filename, safe = TRUE, verbose = TRUE,
                            ext = '.rds'){
    if(!is.null(filename) && safe){ #save file
        filename <- gsub('.rds', "", filename)
        filename0 <- filename
        count <- 1L
        # create a new file name if old one exists, and throw warning
        while(TRUE){
            filename <- paste0(filename, ext)
            if(file.exists(filename)){
                filename <- paste0(filename0, '-', count)
                count <- count + 1L
            } else break
        }
        if(count > 1L)
            if(verbose && safe)
                message(paste0('\nWARNING:\n', filename0, 'existed in the working directory.
                               Using a unique file name instead.\n'))
    }
    filename
}

load_packages <- function(packages){
    if(!is.null(packages)){
        for(pack in packages){
            available <- suppressWarnings(require(substitute(pack), character.only=TRUE,
                    quietly=TRUE, warn.conflicts=FALSE))
            if(!available)
                stop(sprintf("Package \'%s\' is not available. Please install.", pack),
                     call.=FALSE)
        }
    }
    invisible(NULL)
}

get_packages <- function(packages){
    sapply(packages, function(x) as.character(packageVersion(x)))
}

# base-code borrowed and modified from pbapply
timeFormater_internal <- function(time, decimals = TRUE){
    dec <- time - floor(time)
    time <- floor(time - dec)
    dec <- round(dec, 2)
    sec <- round(time %% 60)
    if(decimals) sec <- sec + dec
    time <- floor(time / 60)
    minutes <- floor(time %% 60)
    time <- floor(time / 60)
    days <- floor(time / 24)
    time <- floor(time %% 24)
    hours <- floor(time %% 60)
    resTime <- ""
    if (days > 0)
        resTime <- sprintf("%02id ", days)
    if (hours > 0 || days > 0)
        resTime <- paste(resTime, sprintf("%02ih ", hours), sep = "")
    if (minutes > 0 || hours > 0 || days > 0)
        resTime <- paste(resTime, sprintf("%02im ", minutes), sep = "")
    resTime <- if(decimals) paste0(resTime, sprintf("%.2fs", sec))
    else paste0(resTime, sprintf("%02is", sec))
    resTime
}

print_progress <- function(row, trow, stored_time, RAM, progress,
                           condition, replications){
    if(progress) cat('\n')
    tmp <- as.list(subset(condition, select=colnames(condition) != "ID"))
    nms <- names(tmp)
    nms2 <- do.call(c, lapply(tmp, as.character))
    wdth <- 85 - 13
    condstring <- paste0(nms, '=', nms2, collapse=', ')
    if(nchar(condstring) > wdth){
        nms <- abbreviate(nms, minlength = 6)
        condstring <- paste0(nms, '=', nms2, collapse=', ')
        if(nchar(condstring) > wdth){
            nms2 <- abbreviate(nms2)
            condstring <- paste0(nms, '=', nms2, collapse=', ')
        }
    }
    if(RAM != "") RAM <- sprintf(';   RAM Used: %s;', RAM)
    if(row == 1 && trow == 1)
        cat(sprintf('\rReplications: %i%s   ', replications, RAM))
    else
        cat(sprintf('\rDesign: %i/%i;   Replications: %i%s   Total Time: %s ',
                    row, trow, replications, RAM, timeFormater_internal(sum(stored_time))))
    cat(sprintf('\n Conditions: %s\n', condstring))
    if(progress) cat('\r')
    utils::flush.console()
    invisible(NULL)
}

myundebug <- function(fun) if(isdebugged(fun)) undebug(fun)

#' Suppress verbose function messages
#'
#' This function is used to suppress information printed from external functions
#' that make internal use of \code{\link{message}} and \code{\link{cat}}, which
#' provide information in interactive R sessions. For simulations, the session
#' is not interactive, and therefore this type of output should be suppressed.
#' For similar behaviour for suppressing warning messages, see
#' \code{\link{manageWarnings}}.
#'
#' @param ... the functional expression to be evaluated
#'
#' @param cat logical; also capture calls from \code{\link{cat}}? If
#'   \code{FALSE} only \code{\link{message}} will be suppressed
#'
#' @param keep logical; return a character vector of the messages/concatenate
#'   and print strings as an attribute to the resulting object from \code{expr(...)}?
#'
#' @param attr.name attribute name to use when \code{keep = TRUE}
#'
#' @seealso \code{\link{manageWarnings}}
#'
#' @export
#'
#' @references
#'
#' Chalmers, R. P., & Adkins, M. C.  (2020). Writing Effective and Reliable Monte Carlo Simulations
#' with the SimDesign Package. \code{The Quantitative Methods for Psychology, 16}(4), 248-280.
#' \doi{10.20982/tqmp.16.4.p248}
#'
#' Sigal, M. J., & Chalmers, R. P. (2016). Play it again: Teaching statistics with Monte
#' Carlo simulation. \code{Journal of Statistics Education, 24}(3), 136-156.
#' \doi{10.1080/10691898.2016.1246953}
#'
#' @examples
#'
#' myfun <- function(x, warn=FALSE){
#'    message('This function is rather chatty')
#'    cat("It even prints in different output forms!\n")
#'    message('And even at different....')
#'    cat("...times!\n")
#'    if(warn)
#'      warning('It may even throw warnings!')
#'    x
#' }
#'
#' out <- myfun(1)
#' out
#'
#' # tell the function to shhhh
#' out <- quiet(myfun(1))
#' out
#'
#' # which messages are suppressed? Extract stored attribute
#' out <- quiet(myfun(1), keep = TRUE)
#' attr(out, 'quiet.messages')
#'
#' # Warning messages still get through (see manageWarnings(suppress)
#' #  for better alternative than using suppressWarnings())
#' out2 <- myfun(2, warn=TRUE) |> quiet() # warning gets through
#' out2
#'
#' # suppress warning message explicitly, allowing others to be raised if present
#' myfun(2, warn=TRUE) |> quiet() |>
#'    manageWarnings(suppress='It may even throw warnings!')
#'
quiet <- function(..., cat=TRUE, keep=FALSE, attr.name='quiet.messages'){
    fun <- function(x) eval(x)
    capts <- NULL
    mess <- if(cat)
        testthat::capture_messages(
            testthat::capture_output_lines(ret <- fun(...)) -> capts)
    else testthat::capture_messages(ret <- fun(...))
    if(keep)
        attr(ret, attr.name) <- c(message.=mess, cat.=capts)
    ret
}

#' Auto-named Concatenation of Vector or List
#'
#' This is a wrapper to the function \code{\link{c}}, however names the respective elements
#' according to their input object name. For this reason, nesting \code{nc()} calls
#' is not recommended (joining independent \code{nc()} calls via \code{c()}
#' is however reasonable).
#'
#' @param ... objects to be concatenated
#'
#' @param use.names logical indicating if \code{names} should be preserved (unlike \code{\link{c}},
#'   default is \code{FALSE})
#'
#' @param error.on.duplicate logical; if the same object name appears in the returning object
#'   should an error be thrown? Default is \code{TRUE}
#'
#' @export
#'
#' @references
#'
#' Chalmers, R. P., & Adkins, M. C.  (2020). Writing Effective and Reliable Monte Carlo Simulations
#' with the SimDesign Package. \code{The Quantitative Methods for Psychology, 16}(4), 248-280.
#' \doi{10.20982/tqmp.16.4.p248}
#'
#' Sigal, M. J., & Chalmers, R. P. (2016). Play it again: Teaching statistics with Monte
#' Carlo simulation. \code{Journal of Statistics Education, 24}(3), 136-156.
#' \doi{10.1080/10691898.2016.1246953}
#'
#' @examples
#'
#' A <- 1
#' B <- 2
#' C <- 3
#'
#' names(C) <- 'LetterC'
#'
#' # compare the following
#' c(A, B, C) # unnamed
#'
#' nc(A, B, C) # named
#' nc(this=A, B, C) # respects override named (same as c() )
#' nc(this=A, B, C, use.names = TRUE) # preserve original name
#'
#' \dontrun{
#' # throws errors if names not unique
#' nc(this=A, this=B, C)
#' nc(LetterC=A, B, C, use.names=TRUE)
#' }
#'
#' # poor input choice names
#' nc(t.test(c(1:2))$p.value, t.test(c(3:4))$p.value)
#'
#' # better to explicitly provide name
#' nc(T1 = t.test(c(1:2))$p.value,
#'    T2 = t.test(c(3:4))$p.value)
#'
#' # vector of unnamed inputs
#' A <- c(5,4,3,2,1)
#' B <- c(100, 200)
#'
#' nc(A, B, C) # A's and B's numbered uniquely
#' c(A, B, C)  # compare
#' nc(beta=A, B, C) # replacement of object name
#'
#' # retain names attributes (but append object name, when appropriate)
#' names(A) <- letters[1:5]
#' nc(A, B, C)
#' nc(beta=A, B, C)
#' nc(A, B, C, use.names=TRUE)
#'
#' # mix and match if some named elements work while others do not
#' c( nc(A, B, use.names=TRUE), nc(C))
#'
#' \dontrun{
#' # error, 'b' appears twice
#' names(B) <- c('b', 'b2')
#' nc(A, B, C, use.names=TRUE)
#' }
#'
#' # List input
#' A <- list(1)
#' B <- list(2:3)
#' C <- list('C')
#'
#' names(C) <- 'LetterC'
#'
#' # compare the following
#' c(A, B, C) # unnamed
#'
#' nc(A, B, C) # named
#' nc(this=A, B, C) # respects override named (same as c() and list() )
#' nc(this=A, B, C, use.names = TRUE) # preserve original name
#'
#'
nc <- function(..., use.names=FALSE, error.on.duplicate = TRUE){
    dots <- list(...)
    len <- sapply(dots, length)
    object <- as.list(substitute(list(...)))[-1L]
    nms <- sapply(object, function(x) paste0(as.character(x), collapse='_'))
    nms[names(nms) != ""] <- names(nms[names(nms) != ""])
    if(any(len > 1L)){
        nms <- as.list(nms)
        for(i in length(nms):1L){
            if(len[i] > 1L)
                nms[[i]] <- paste0(rep(nms[[i]], len[i]),
                                   if(!is.null(names(dots[[i]]))) "." else NULL,
                                if(is.null(names(dots[[i]]))) 1L:len[i]
                                else names(dots[[i]]))
        }
        nms <- do.call(c, nms)
    }
    if(use.names){
        tmp <- do.call(c, lapply(dots, function(x){
            ret <- names(x)
            if(is.null(ret)) ret <- rep(NA, length(x))
            ret
        }))
        nms[!is.na(tmp)] <- tmp[!is.na(tmp)]
    }
    nms <- gsub("\\$\\_", "", nms)
    if(error.on.duplicate)
        if(any(duplicated(nms)))
            stop(sprintf('Vector/list contains the following duplicated names: %s',
                         paste0(nms[duplicated(nms)], collapse=', ')),
                 call.=FALSE)
    ret <- c(...)
    names(ret) <- nms
    ret
}

#' Check whether package versions are as expected
#'
#' Check if a specific version of a package is installed, and throw an error if not what was
#' anticipated. Particularly useful when submitting jobs to cluster compute nodes where
#' package versions may be been updated inadvertently.
#'
#' @param ... character vectors indicating package and version to check,
#'   of the form \code{"package OPERATOR version"}, where
#'   \code{"OPERATOR"} is one of R's logical expressions (e.g., \code{"dplyr == 1.2.1"}).
#'   Can contain one or more expressions to evaluate
#'
#' @return invisible return of the package names as a character vector
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#'   CheckPackages('dplyr == 1.2.1') # fails if not exact version
#'   CheckPackages('dplyr <= 1.2.1') # fails if not equal to or less than
#'   CheckPackages('dplyr >= 1.2.1') # allows specific version or higher
#'
#'   CheckPackages('dplyr <= 1.2.1', 'mirt >= 1.45.1')  # multiple checks
#'
#' }
CheckPackages <- function(...){
    packageVersion <- list(...)
    pack.out <- sapply(packageVersion, \(pack){
        pack <- gsub(" ", "", pack, fixed = TRUE)
        if(grepl('==', pack, fixed=TRUE)){
            split <- strsplit(pack, "==", fixed = TRUE)[[1]]
            if(!package_version(packageVersion(split[1])) == split[2])
                stop(sprintf("%s version is not equal to %s", split[1], split[2]),
                     call.=FALSE)
        } else if(grepl('>', pack, fixed=TRUE)){
            if(grepl('>=', pack, fixed=TRUE)){
                split <- strsplit(pack, ">=", fixed = TRUE)[[1]]
                if(!package_version(packageVersion(split[1])) >= split[2])
                    stop(sprintf("%s version is not greater than or equal to %s", split[1], split[2]),
                         call.=FALSE)
            } else {
                split <- strsplit(pack, ">", fixed = TRUE)[[1]]
                if(!package_version(packageVersion(split[1])) > split[2])
                    stop(sprintf("%s version is not greater than %s", split[1], split[2]), call.=FALSE)
            }
        } else if(grepl('<', pack, fixed=TRUE)){
            if(grepl('<=', pack, fixed=TRUE)){
                split <- strsplit(pack, "<=", fixed = TRUE)[[1]]
                if(!package_version(packageVersion(split[1])) <= split[2])
                    stop(sprintf("%s version is not less than or equal to %s", split[1], split[2]),
                         call.=FALSE)
            } else {
                split <- strsplit(pack, "<", fixed = TRUE)[[1]]
                if(!package_version(packageVersion(split[1])) < split[2])
                    stop(sprintf("%s version is not less than %s", split[1], split[2]), call.=FALSE)
            }
        } else stop('package specification is malformed')
        split[1]
    })
    invisible(pack.out)
}

isList <- function(x) !is.data.frame(x) && is.list(x)

reduceTable <- function(tab){
    tab <- dplyr::bind_rows(tab)
    uniq <- sort(unique(tab$x))
    reps <- val <- numeric(length(uniq))
    for(i in seq_len(length(val))){
        tmp <- tab[uniq[i] == tab$x, , drop=FALSE]
        reps[i] <- sum(tmp$reps)
        val[i] <- sum(as.numeric(tmp$y) * as.numeric(tmp$reps) / reps[i])
    }
    reduced <- data.frame(y=val, x=uniq, reps=reps)
    reduced

}

sim_results_check <- function(sim_results, return_list = FALSE){
    if(is(sim_results, 'try-error')){
        warning(c("Summarise() should not throw errors; please inspect. \nError raised:\n    ",
                  sim_results), call.=FALSE)
        return(c(SUMMARISE_ERROR=NA))
    }
    if(is.data.frame(sim_results) || is.matrix(sim_results)){
        if(nrow(sim_results) > 1L){
            sim_results <- list(sim_results)
        } else {
            nms <- colnames(sim_results)
            sim_results <- as.numeric(sim_results)
            names(sim_results) <- nms
        }
    }
    if(isList(sim_results)){
        if(length(sim_results) > 1L && is.null(names(sim_results)))
            stop("List elements must be named in Summarise() definition",
                 call.=FALSE)
        if(return_list) return(sim_results)
        ret <- numeric(0)
        attr(ret, 'summarise_list') <- sim_results
        return(ret)
    }
    if(length(sim_results) == 1L){
        if(is.null(names(sim_results)))
            names(sim_results) <- 'value'
        if(!is.vector(sim_results) || is.null(names(sim_results)))
            stop('summarise() must return a named vector or data.frame object with 1 row',
                 call.=FALSE)
    }
    sim_results
}

unwind_apply_wind.list <- function(lst, mat, fun, ...){
    long_list <- do.call(rbind, lapply(lst, as.numeric))
    long_mat <- if(!is.null(mat)) as.numeric(mat) else NULL
    ret <- fun(long_list, long_mat, ...)
    if(!is.null(mat)){
        was_matrix <- is.matrix(mat)
        if(was_matrix){
            ret <- matrix(ret, nrow(mat), ncol(mat))
            rownames(ret) <- rownames(mat)
            colnames(ret) <- colnames(mat)
        } else names(ret) <- names(mat)
    }
    ret
}

lapply_timer <- function(X, FUN, max_time, max_RAM, ...){
    if(is.finite(max_time)){
        ret <- vector('list', length(X))
        total <- max_time
        elapsed <- 0
        time_left <- total
        for(i in 1L:length(ret)){
            st <- proc.time()['elapsed']
            val <- R.utils::withTimeout(FUN(i, ...),
                                        timeout = time_left,
                                        onTimeout = 'warning')
            elapsed <- elapsed + proc.time()['elapsed'] - st
            time_left <- total - elapsed
            ret[[i]] <- val
            if(time_left <= 0){
                message(sprintf(c("Simulation terminated due to max_time constraint",
                                " (%i/%i replications evaluated)."), i, length(ret)))
                ret <- ret[1L:i]
                break
            }
            if(is.finite(max_RAM) && object.size(ret) > max_RAM){
                message(sprintf(c("Simulation terminated due to max_RAM constraint",
                                  " (%i/%i replications evaluated)."), i, length(ret)))
                ret <- ret[1L:i]
                break
            }
        }
    } else {
        ret <- vector('list', length(X))
        for(i in 1L:length(ret)){
            val <- FUN(i, ...)
            ret[[i]] <- val
            if(is.finite(max_RAM) && object.size(ret) > max_RAM){
                message(sprintf(c("Simulation terminated due to max_RAM constraint",
                                  " (%i/%i replications evaluated)."), i, length(ret)))
                ret <- ret[1L:i]
                break
            }
        }
    }
    ret
}

combined_Analyses <- function(condition, dat, fixed_objects){
    if(!is.null(.SIMDENV$ANALYSE_FUNCTIONS)){
        ANALYSE_FUNCTIONS <- .SIMDENV$ANALYSE_FUNCTIONS
        TRY_ALL_ANALYSE <- .SIMDENV$TRY_ALL_ANALYSE
    }
    nfuns <- length(ANALYSE_FUNCTIONS)
    ret <- vector('list', nfuns)
    nms <- names(ANALYSE_FUNCTIONS)
    names(ret) <- nms
    if(is.null(nms)) nms <- 1L:nfuns
    for(i in nms){
        tried <- try(ANALYSE_FUNCTIONS[[i]](condition=condition, dat=dat,
                                            fixed_objects=fixed_objects), silent=TRUE)
        if(is(tried, 'try-error')){
            if(tried == 'Error : ANALYSEIF RAISED ERROR\n')
                tried <- NULL
            else if(!TRY_ALL_ANALYSE) return(tried)
        }
        ret[[i]] <- tried
    }
    if(TRY_ALL_ANALYSE){
        try_error <- sapply(ret, function(x) is(x, 'try-error'))
        if(any(try_error)){
            msg <- paste0(names(ANALYSE_FUNCTIONS)[try_error], ".ERROR:   ", as.character(ret[try_error]))
            if(length(msg) > 1L)
                msg[1L] <- sprintf("%i INDEPENDENT ERRORS THROWN:   %s", length(msg), msg[1L])
            msg <- paste0(msg, collapse = '  ')
            ret <- try(stop(msg), silent = TRUE)
            ret <- gsub("Error in try\\(stop\\(msg\\), silent = TRUE\\) : \\\n  ", "", ret)
            return(ret)
        }
    }
    if(all(sapply(ret, function(x) is.numeric(x) ||
                  (is.data.frame(x) && nrow(x) == 1L))))
        ret <- unlist(ret)
    ret
}

combined_Generate <- function(condition, fixed_objects){
    if(!is.null(.SIMDENV$GENERATE_FUNCTIONS))
        GENERATE_FUNCTIONS <- .SIMDENV$GENERATE_FUNCTIONS
    nfuns <- length(GENERATE_FUNCTIONS)
    ret <- vector('list', nfuns)
    nms <- names(GENERATE_FUNCTIONS)
    names(ret) <- nms
    if(is.null(nms)) nms <- 1L:nfuns
    for(i in nms){
        tried <- try(GENERATE_FUNCTIONS[[i]](condition=condition,
                                            fixed_objects=fixed_objects), silent=TRUE)
        if(is(tried, 'try-error')){
            if(tried == 'Error : GENERATEIF RAISED ERROR\n')
                tried <- NULL
        } else ret <- tried
    }
    if(is.null(ret))
        stop('No data was generated for supplied condition. Please fix', call.=FALSE)
    ret
}

toTabledResults <- function(results){
    tabled_results <- if(is.data.frame(results[[1]]) && nrow(results[[1L]]) == 1L){
        dplyr::bind_rows(results)
    } else if((is.data.frame(results[[1]]) && nrow(results[[1]]) > 1L) || is.list(results[[1L]])){
        results
    } else {
        as.matrix(dplyr::bind_rows(as.data.frame(do.call(rbind, results))))
    }
    tabled_results
}

stackResults <- function(results){
    if(!is.list(results[[1L]]) || (is.data.frame(results[[1L]]) &&
                                   nrow(results[[1L]]) == 1L)){
        old_nms <- names(results[[1L]])
        results <- as.data.frame(do.call(rbind, results))
        if(length(unique(colnames(results))) != ncol(results) && ncol(results) > 1L)
            stop('Object of results returned from analyse must have unique names', call.=FALSE)
        rownames(results) <- NULL
        if(ncol(results) == 1L && is.null(old_nms)) results <- results[,1]
    }
    results
}

SimSolveData <- function(burnin, full = TRUE){
    pick <- !sapply(.SIMDENV$stored_results, is.null)
    pick[1L:burnin] <- FALSE
    if(!any(pick))
        return(data.frame(y=numeric(0), IV=numeric(0), weights=numeric(0)))
    if(full){
        DV <- do.call(c, .SIMDENV$stored_results[pick])
        IV <- rep(.SIMDENV$stored_medhistory[pick],
                  times=sapply(.SIMDENV$stored_results[pick], length))
        ret <- data.frame(y=DV, x=IV, weights=1)
    } else {
        ret <- do.call(rbind, .SIMDENV$stored_history[pick])
        ret$weights <- 1/sqrt(ret$reps)
    }
    ret
}

SimSolveUniroot <- function(SimMod, b, interval, max.interval, median, CI=NULL){
    f.root <- function(x, b)
        predict(SimMod, newdata = data.frame(x=x), type = 'response') - b
    res <- try(uniroot(f.root, b=b, interval = interval), silent = TRUE)
    if(is(res, 'try-error')){
        org.interval <- interval
        # in case original interval is poor for interpolation
        interval <- max.interval
        for(i in seq_len(20L)){
            if(grepl('end points not of opposite sign', res)){
                diff <- abs(interval - median)
                interval[which.max(diff)] <- mean(c(median, interval[which.max(diff)]))
                res <- try(uniroot(f.root, b=b, interval = interval), silent = TRUE)
                if(!is(res, 'try-error')) break
            }
        }
    }
    if(is(res, 'try-error')) return(c(NA, NA, NA))
    root <- res$root
    abias <- bias(root, median, type = 'abs_relative')
    if(abias > .5) root <- median
    ci <- c(NA, NA)
    if(!is.null(CI)){
        preds <- predict(SimMod, newdata = data.frame(x=root),
                      se.fit=TRUE, type = 'link')
        ci <- SimMod$family$linkinv(preds$fit + qnorm(CI) * preds$se.fit)
    }
    c(root, ci)
}

collect_unique <- function(x){
    if(any(duplicated(colnames(x)))){
        uniq <- unique(colnames(x))
        for(u in uniq){
            pick <- colnames(x) %in% u
            if(sum(pick) == 1L) next
            whc <- sort(which(pick))
            tmp <- rowSums(x[,pick, drop=FALSE], na.rm = TRUE)
            tmp <- ifelse(tmp == 0, NA, tmp)
            x[[whc[1L]]] <- tmp
            x[whc[2L:length(whc)]] <- NULL
        }
    }
    x
}

bisection <- function (f, interval, ..., tol = 0.001, maxiter = 100,
                       f.lower = NULL, f.upper = NULL, check = FALSE)
{
    lower <- interval[1L]
    upper <- interval[2L]
    iter <- 0L
    if(check){
        if(is.null(f.lower)) f.lower <- f(lower, ...)
        if(is.null(f.upper)) f.upper <- f(upper, ...)
        stopifnot("No root in specified interval" = f.lower * f.upper < 0)
    } else {
        if(is.null(f.lower)) f.lower <- -Inf
        if(is.null(f.upper)) f.upper <- Inf
    }
    if(f.lower > f.upper){
        tmp <- lower
        lower <- upper
        upper <- tmp
        tmp <- f.lower
        f.lower <- f.upper
        f.upper <- tmp
    }
    false_converge <- FALSE
    for(i in 1L:maxiter){
        iter <- iter + 1L
        mid <- (lower + upper)/2
        f.mid <- f(mid, ...)
        if(f.mid < f.lower || f.mid > f.upper){
            false_converge <- TRUE
            break
        }
        if (isTRUE(f.lower * f.mid > 0)){
            lower <- mid
            f.lower <- f.mid
        } else {
            upper <- mid
            f.upper <- f.mid
        }
        if(abs(lower - upper) < tol) break
    }
    root <- (lower + upper)/2
    list(root=root, f.root=f(root, ...), iter=i,
         terminated_early=i < maxiter,
         false_converge=false_converge)
}

RAM_used <- function(format=TRUE){
    # borrowed and modified from pryr::node_size(), 13-06-2023
    bit <- 8L * .Machine$sizeof.pointer
    if (!(bit == 32L || bit == 64L)) {
        stop("Unknown architecture", call. = FALSE)
    }
    val <- if (bit == 32L) 28L else 56L
    # end borrowed portion
    bytes <- sum(gc()[, 1] * c(val, 8))
    size <- structure(bytes, class="object_size")
    if(!format) return(size)
    format(size, 'MB')
}

clip_names <- function(vec, maxchar = 150L){
    names(vec) <- strtrim(names(vec), width=maxchar)
    vec
}

#' Form Column Standard Deviation and Variances
#'
#' Form column standard deviation and variances for numeric arrays (or data frames).
#'
#' @param x an array of two dimensions containing numeric, complex, integer or logical values,
#'   or a numeric data frame
#'
#' @param na.rm logical; remove missing values in each respective column?
#'
#' @param unname logical; apply \code{\link{unname}} to the results to remove any variable
#'   names?
#'
#' @seealso \code{\link{colMeans}}
#'
#' @export
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#'
#' @examples
#'
#' results <- matrix(rnorm(100), ncol=4)
#' colnames(results) <- paste0('stat', 1:4)
#'
#' colVars(results)
#' colSDs(results)
#'
#' results[1,1] <- NA
#' colSDs(results)
#' colSDs(results, na.rm=TRUE)
#' colSDs(results, na.rm=TRUE, unname=TRUE)
#'
colVars <- function(x, na.rm=FALSE, unname=FALSE){
    ret <- apply(x, 2L, FUN = var, na.rm=na.rm)
    if(unname) ret <- unname(ret)
    ret
}

#' @export
#' @rdname colVars
colSDs <- function(x, na.rm=FALSE, unname=FALSE){
    sqrt(colVars(x=x, na.rm=na.rm, unname=unname))
}

pickReps <- function(replications, iter){
    ret <- if(iter > length(replications))
        max(replications) else replications[iter]
    ret
}

set_seed <- function(seed){
    if(is.list(seed)) .GlobalEnv$.Random.seed <- seed[[1L]]
    else set.seed(seed)
    invisible(NULL)
}

recvResult_fun <- utils::getFromNamespace("recvResult", "parallel")

#' Set RNG sub-stream for  Pierre L'Ecuyer's RngStreams
#'
#' Sets the sub-stream RNG state within for Pierre L'Ecuyer's (1999)
#' algorithm. Should be used within distributed array jobs
#' after suitable L'Ecuyer's (1999) have been distributed to each array, and
#' each array is further defined to use multi-core processing. See
#' \code{\link[parallel]{clusterSetRNGStream}} for further information.
#'
#' @param seed An integer vector of length 7 as given by \code{.Random.seed} when
#'   the L'Ecuyer-CMR RNG is in use. See\code{\link{RNG}} for the valid values
#' @param cl A cluster from the \code{parallel} package, or
#'   (if \code{NULL}) the registered cluster
#' @return invisible NULL
#' @export
#'
#'
clusterSetRNGSubStream <- function(cl, seed){
    nc <- length(cl)
    seeds <- vector("list", nc)
    seeds[[1L]] <- seed[[1L]]
    for (i in seq_len(nc - 1L)) seeds[[i + 1L]] <-
        parallel::nextRNGSubStream(seeds[[i]])
    for (i in seq_along(cl)) {
        expr <- substitute(assign(".Random.seed", seed, envir = .GlobalEnv),
                           list(seed = seeds[[i]]))
        sendCall.imp(cl[[i]], eval, list(expr))
    }
    checkForRemoteErrors.imp(lapply(cl, recvResult_fun))
    invisible()
}

sendCall.imp <- utils::getFromNamespace('sendCall', 'parallel')
checkForRemoteErrors.imp <- utils::getFromNamespace('checkForRemoteErrors',
                                                    'parallel')
recvOneResult.imp <- utils::getFromNamespace('recvOneResult', 'parallel')

# Generate one L'Ecuyer-CMRG .Random.seed state per replication. Scalar seeds
# are expanded into successive streams (analogous to clusterSetRNGStream(),
# but over replications rather than nodes), while list seeds (which already
# contain a stream state; see genSeeds(iseed)) are expanded into sub-streams
# within that stream (analogous to clusterSetRNGSubStream()) so that streams
# distributed across conditions/arrays remain non-overlapping
gen_replication_RNGseeds <- function(seed, replications){
    seeds <- vector("list", replications)
    if(is.list(seed)){
        seeds[[1L]] <- seed[[1L]]
        for(i in seq_len(replications - 1L))
            seeds[[i + 1L]] <- parallel::nextRNGSubStream(seeds[[i]])
    } else {
        oldseed <- if(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
        rngkind <- RNGkind()
        RNGkind("L'Ecuyer-CMRG")
        set.seed(seed)
        seeds[[1L]] <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        for(i in seq_len(replications - 1L))
            seeds[[i + 1L]] <- parallel::nextRNGStream(seeds[[i]])
        if(!is.null(oldseed)){
            assign(".Random.seed", oldseed, envir = .GlobalEnv)
        } else {
            RNGkind(rngkind[1L])
            set.seed(NULL)
        }
    }
    seeds
}

# Worker-side receiver for LB_dispatch_job(): the invariant arguments are
# stored once per node (see dynamicClusterLapply()) so that each dispatched
# job only carries its replication index and RNG state
LB_store_args <- function(args){
    .SIMDENV$LB_args <- args
    invisible(NULL)
}

LB_dispatch_job <- function(x){
    if(!is.null(x$rng_seed))
        assign(".Random.seed", x$rng_seed, envir = .GlobalEnv)
    do.call(.SIMDENV$LB_args$fun, c(list(x$index), .SIMDENV$LB_args$args))
}

# Dispatcher-backed mirai daemons for load-balanced dispatching. The
# miraiCluster objects from mirai::make_cluster() are hard-coded to
# dispatcher=FALSE, in which case tasks are distributed round-robin at the
# socket level and can queue behind a busy daemon even when other daemons are
# free, making true load balancing impossible through the parallel-API
# compatibility layer. mirai's dispatcher, by contrast, holds queued tasks
# centrally and assigns each to the next free daemon
make_mirai_dispatcher <- function(ncores){
    # unique per-call profile so a re-entrant runSimulation() (e.g., launched
    # from a prepare()/summarise() definition) cannot tear down the daemons
    # belonging to a still-active outer run
    count <- if(is.null(.SIMDENV$LB_profile_count)) 1L
        else .SIMDENV$LB_profile_count + 1L
    .SIMDENV$LB_profile_count <- count
    profile <- sprintf('SimDesign.LB.%i', count)
    # clear any stale daemons left over from an interrupted previous run
    try(mirai::daemons(0L, .compute = profile), silent = TRUE)
    mirai::daemons(ncores, dispatcher = TRUE, .compute = profile)
    structure(list(profile = profile, ncores = ncores),
              class = 'SimDesignMiraiLB')
}

stop_mirai_dispatcher <- function(cl)
    try(mirai::daemons(0L, .compute = cl$profile), silent = TRUE)

length.SimDesignMiraiLB <- function(x) x$ncores

# everywhere()-based analogue of parallel::clusterExport(); broadcast is
# guaranteed to reach every daemon, which per-node sends are not once the
# dispatcher (rather than the node object) decides task placement
mirai_export <- function(cl, varlist, envir){
    if(!length(varlist)) return(invisible(NULL))
    objs <- lapply(varlist, function(nm) get(nm, envir = envir))
    names(objs) <- varlist
    # everywhere() assigns its ... objects into the daemon global environments
    # (persistently), which is also where the user's exports live; use a
    # dotted carrier name that neither ls() discovery nor user code will
    # collide with, and remove it in the same evaluation
    mirai::everywhere({
        list2env(.sd_exports, envir = globalenv())
        rm('.sd_exports', envir = globalenv())
    }, .sd_exports = objs, .compute = cl$profile)
    invisible(NULL)
}

# native-mirai analogue of dynamicClusterLapply(); all jobs are queued with
# the dispatcher up-front, which streams them to daemons as they free up
mirai_dispatchLapply <- function(cl, X, fun, ..., rng_seeds = NULL, progress = FALSE){
    n <- length(X)
    val <- vector("list", n)
    if(n > 0L){
        # rm() in-expression so the daemon global environments neither clobber
        # same-named user exports nor retain a duplicate copy of the arguments;
        # the working copy lives in each daemon's .SIMDENV until cleared below
        mirai::everywhere({
            .sd_store(.sd_args)
            rm('.sd_store', '.sd_args', envir = globalenv())
        }, .sd_store = LB_store_args, .sd_args = list(fun=fun, args=list(...)),
           .compute = cl$profile)
        jobs <- lapply(seq_len(n), function(job)
            mirai::mirai(RUN(x), RUN = LB_dispatch_job,
                         x = list(index=X[[job]],
                                  rng_seed=if(is.null(rng_seeds)) NULL
                                           else rng_seeds[[job]]),
                         .compute = cl$profile))
        pb <- NULL
        if(progress){
            pb <- pbapply::startpb(0, n)
            on.exit(pbapply::closepb(pb), add = TRUE)
        }
        # Wait with a daemon-liveness watchdog rather than blocking in
        # call_mirai(). The dispatcher resolves *in-flight* tasks to error
        # values when their daemon dies, but tasks still centrally queued
        # when the last daemon dies would remain unresolved forever (e.g.,
        # replications that crash their worker process through segfaulting
        # compiled code or OOM kills), hanging the master R session. The
        # zero-connection state must persist across consecutive checks
        # before aborting since daemons may not have connected yet when
        # the first jobs are queued
        zero_count <- 0L
        repeat{
            unres <- vapply(jobs, mirai::unresolved, logical(1L))
            if(progress) pbapply::setpb(pb, n - sum(unres))
            if(!any(unres)) break
            # a resolved try-error is fatal (mainsim() internally retries the
            # recoverable errors), so cancel the still-queued jobs rather than
            # letting every remaining replication fail through max_errors
            # generate/analyse attempts of its own before the failure surfaces
            if(any(vapply(jobs[!unres], function(m) inherits(m$data, 'try-error'),
                          logical(1L)))){
                for(m in jobs[unres]) try(mirai::stop_mirai(m), silent = TRUE)
                break
            }
            connections <- try(mirai::status(.compute = cl$profile)$connections,
                               silent = TRUE)
            zero_count <- if(!is(connections, 'try-error') && isTRUE(connections == 0L))
                zero_count + 1L else 0L
            if(zero_count >= 100L){    # ~10 seconds with no live daemons
                for(m in jobs[unres]) try(mirai::stop_mirai(m), silent = TRUE)
                stop(sprintf('all %i mirai daemons terminated while %i replication(s) remained unevaluated',
                             cl$ncores, sum(unres)), call.=FALSE)
            }
            Sys.sleep(.1)
        }
        val <- lapply(jobs, function(m) mirai::call_mirai(m)$data)
        mirai::everywhere({
            .sd_store(NULL)
            rm('.sd_store', envir = globalenv())
        }, .sd_store = LB_store_args, .compute = cl$profile)
        # report the original failure rather than the errorValue of a job that
        # was cancelled above or lost with a dying daemon
        real_errs <- vapply(val, inherits, logical(1L), what = 'miraiError')
        if(any(real_errs))
            stop('one node produced an error: ',
                 as.character(val[[which(real_errs)[1L]]]), domain = NA, call. = FALSE)
        # bare errorValues mean the daemon died mid-evaluation; phrase the
        # message so the fatal-termination handler can extract the reason
        lost <- vapply(val, function(x) inherits(x, 'errorValue') &&
                           !inherits(x, 'miraiError'), logical(1L))
        if(any(lost))
            stop(paste0('One or more parallel daemons terminated while evaluating ',
                        'replications. \n\nLast error message was: \n\n  ',
                        'daemon connection lost (errorValue ',
                        as.integer(val[[which(lost)[1L]]]),
                        '), e.g., a worker crash via segfault or out-of-memory kill'),
                 call. = FALSE)
    }
    checkForRemoteErrors.imp(val)
}

# Load-balanced alternative to parallel::parLapply() with optional
# per-replication RNG states and pbapply-based progress reporting. Each
# element of X is dispatched to a node as soon as that node finishes its
# previous job (cf. parallel:::dynamicClusterApply()), so replications with
# heterogeneous run-times no longer leave the remaining workers idle, which
# occurs with statically scheduled dispatches (parLapply()) and
# batch-synchronized progress dispatches (pblapply(cl=))
dynamicClusterLapply <- function(cl, X, fun, ..., rng_seeds = NULL, progress = FALSE){
    n <- length(X)
    p <- length(cl)
    val <- vector("list", n)
    if(n > 0L && p > 0L){
        parallel::clusterCall(cl, LB_store_args, list(fun=fun, args=list(...)))
        # Only clean up the worker-side state after a successful drain. If a
        # worker died mid-run the sockets are in an undefined state, and the
        # write performed here would consume the single post-close write that
        # succeeds silently, poisoning the later stopCluster() call with an
        # uncaught SIGPIPE error that discards runSimulation()'s fully
        # assembled return value
        drained <- FALSE
        on.exit(if(drained)
            try(parallel::clusterCall(cl, LB_store_args, NULL), silent = TRUE),
            add = TRUE)
        submit <- function(node, job)
            sendCall.imp(cl[[node]], LB_dispatch_job,
                         list(list(index=X[[job]],
                                   rng_seed=if(is.null(rng_seeds)) NULL
                                            else rng_seeds[[job]])),
                         tag = job)
        pb <- NULL
        if(progress){
            pb <- pbapply::startpb(0, n)
            on.exit(pbapply::closepb(pb), add = TRUE)
        }
        outstanding <- min(n, p)
        for(i in seq_len(outstanding)) submit(i, i)
        nxt <- outstanding + 1L
        done <- 0L
        had_error <- FALSE
        while(outstanding > 0L){
            d <- recvOneResult.imp(cl)
            outstanding <- outstanding - 1L
            done <- done + 1L
            val[d$tag] <- list(d$value)
            # a try-error here is fatal (mainsim() internally retries recoverable
            # errors), so stop dispatching new jobs and drain the active ones
            if(inherits(d$value, 'try-error')) had_error <- TRUE
            if(!had_error && nxt <= n){
                submit(d$node, nxt)
                nxt <- nxt + 1L
                outstanding <- outstanding + 1L
            }
            if(progress) pbapply::setpb(pb, done)
        }
        drained <- TRUE
    }
    checkForRemoteErrors.imp(val)
}

valid_results <- function(x)
    is(x, 'numeric') || is(x, 'data.frame') || is(x, 'list') || is(x, 'logical') || is(x, 'try-error')

#' Generate random seeds
#'
#' Generate seeds to be passed to \code{runSimulation}'s \code{seed} input. Values
#' are sampled from 1 to 2147483647, or are generated using L'Ecuyer-CMRG's (2002)
#' method (returning either a list if \code{arrayID} is omitted, or the specific
#' row value from this list if \code{arrayID} is included).
#'
#' @param design design matrix that requires a unique seed per condition, or
#'   a number indicating the number of seeds to generate. Default generates one
#'   number
#'
#' @param iseed the initial \code{set.seed} number used to generate a sequence
#'   of independent seeds according to the L'Ecuyer-CMRG (2002) method. This
#'   is recommended whenever quality random number generation is required
#'   across similar (if not identical) simulation jobs
#'   (e.g., see \code{\link{runArraySimulation}}). If \code{arrayID} is not
#'   specified then this will return a list of the associated seed for the
#'   full \code{design}
#'
#' @param arrayID (optional) single integer input corresponding to the specific
#'   row in the \code{design} object when using the \code{iseed} input.
#'   This is used in functions such as \code{\link{runArraySimulation}}
#'   to pull out the specific seed rather than manage a complete list, and
#'   is therefore more memory efficient
#'
#' @param old.seeds (optional) vector or matrix of last seeds used in
#'   previous simulations to avoid repeating the same seed on a subsequent run.
#'   Note that this approach should be used sparingly as seeds set more frequently
#'   are more likely to correlate, and therefore provide less optimal random
#'   number behaviour (e.g., if performing a simulation on two runs to achieve
#'   5000 * 2 = 10,000 replications this is likely reasonable,
#'   but for simulations with 100 * 2 = 200 replications this is more
#'   likely to be sub-optimal).
#'   Length must be equal to the number of rows in \code{design}
#'
#' @export
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#'
#' @examples
#'
#' # generate 1 seed (default)
#' genSeeds()
#'
#' # generate 5 unique seeds
#' genSeeds(5)
#'
#' # generate from nrow(design)
#' design <- createDesign(factorA=c(1,2,3),
#'                        factorB=letters[1:3])
#' seeds <- genSeeds(design)
#' seeds
#'
#' # construct new seeds that are independent from original (use this sparingly)
#' newseeds <- genSeeds(design, old.seeds=seeds)
#' newseeds
#'
#' # can be done in batches too
#' newseeds2 <- genSeeds(design, old.seeds=cbind(seeds, newseeds))
#' cbind(seeds, newseeds, newseeds2) # all unique
#'
#' ############
#' # generate seeds for runArraySimulation()
#' (iseed <- genSeeds())  # initial seed
#' seed_list <- genSeeds(design, iseed=iseed)
#' seed_list
#'
#' # expand number of unique seeds given iseed (e.g., in case more replications
#' # are required at a later date)
#' seed_list_tmp <- genSeeds(nrow(design)*2, iseed=iseed)
#' str(seed_list_tmp) # first 9 seeds identical to seed_list
#'
#' # more usefully for HPC, extract only the seed associated with an arrayID
#' arraySeed.15 <- genSeeds(nrow(design)*2, iseed=iseed, arrayID=15)
#' arraySeed.15
#'
genSeeds <- function(design = 1L, iseed = NULL, arrayID = NULL, old.seeds = NULL){
    if(missing(design)) design <- 1L
    if(is.numeric(design))
        design <- matrix(NA, nrow=design)
    if(is.null(iseed)){
        seed <- rint(nrow(design), min=1L, max = 2147483647L)
        if(!is.null(old.seeds)){
            old.seeds <- as.vector(old.seeds)
            while(TRUE){
                whc <- which(seed %in% old.seeds)
                if(length(whc)){
                    seed[whc] <- rint(nrow(design), min=1L, max = 2147483647L)
                    next
                }
                break
            }
        }
    } else {
        rngkind <- RNGkind()
        RNGkind("L'Ecuyer-CMRG")
        on.exit({RNGkind(rngkind[1L]); set.seed(NULL)})
        seed <- if(!is.null(arrayID)) vector('list', 1L)
            else vector('list', nrow(design))
        set.seed(iseed)
        seed[[1L]] <- .Random.seed
        if(!is.null(arrayID)){
            stopifnot(is.numeric(arrayID) && length(arrayID) == 1L)
            if(arrayID < 1L || arrayID > nrow(design))
                stop('arrayID not associated with valid row in design')
            seed.i <- seed[[1L]]
            if(arrayID > 1L){
                for (i in 2L:arrayID)
                    seed.i <- nextRNGStream(seed.i)
            }
            seed[[1L]] <- seed.i
            attr(seed, 'arrayID') <- arrayID
        } else {
            if(length(seed) > 1L){
                for (i in 2L:length(seed))
                    seed[[i]] <- nextRNGStream(seed[[i - 1L]])
            }
        }
        attr(seed, 'iseed') <- iseed
    }
    seed
}

#' Format time string to suitable numeric output
#'
#' Format time input string into suitable numeric output metric (e.g., seconds).
#' Input follows the \code{SBATCH} utility specifications.
#' Accepted time formats include \code{"minutes"},
#' \code{"minutes:seconds"}, \code{"hours:minutes:seconds"},
#' \code{"days-hours"}, \code{"days-hours:minutes"} and
#' \code{"days-hours:minutes:seconds"}. Alternatively, function can be used to
#' convert numeric input to SBATCH format.
#'
#' For example, \code{time = "60"} indicates a maximum time of 60 minutes,
#' \code{time = "03:00:00"} a maximum time of 3 hours,
#' \code{time = "4-12"} a maximum of 4 days and 12 hours, and
#' \code{time = "2-02:30:00"} a maximum of 2 days, 2 hours and 30 minutes.
#'
#' @param time a character string to be formatted. If a numeric vector is supplied
#' then this will be interpreted as minutes due to character coercion.
#'
#' @param output type of numeric output to convert time into.
#' Currently supported are \code{'sec'} for seconds (default),
#' \code{'min'} for minutes, \code{'hour'}, and \code{'day'}.
#'
#' Alternatively, if \code{time} were numeric then setting \code{output} to
#' \code{'SBATCH'} will return a suitable SBATCH format.
#'
#' @param input if supplied \code{time} is a numeric, indicates what the value
#'   represents. Default assumes the input is in minutes (see \code{output} for
#'   supported values)
#'
#' @param sround function used to round last seconds computation
#'
#' @export
#'
#' @examples
#'
#' # Test cases (outputs in seconds)
#' timeFormater("4-12")        # day-hours
#' timeFormater("4-12:15")     # day-hours:minutes
#' timeFormater("4-12:15:30")  # day-hours:minutes:seconds
#'
#' timeFormater("30")          # minutes
#' timeFormater("30:30")       # minutes:seconds
#' timeFormater("4:30:30")     # hours:minutes:seconds
#'
#' # output in hours
#' timeFormater("4-12", output = 'hour')
#' timeFormater("4-12:15", output = 'hour')
#' timeFormater("4-12:15:30", output = 'hour')
#'
#' timeFormater("30", output = 'hour')
#' timeFormater("30:30", output = 'hour')
#' timeFormater("4:30:30", output = 'hour')
#'
#' # numeric input is understood as minutes
#' timeFormater(42)                            # output in seconds
#' timeFormater(42, input='sec', output='min') # input in sec, output in min
#' timeFormater(42, output='min')              # in-and-output in minutes
#'
#' # convert numeric inputs to SBATCH format
#' timeFormater(60, output='SBATCH')
#' timeFormater(3, output='SBATCH', input='day')
#' timeFormater(7000, output='SBATCH', input='sec')
#' timeFormater(100000, output='SBATCH', input='sec')
#'
#' # rounding seconds
#' timeFormater(1.55555, output='SBATCH', input='sec') # floor default
#' timeFormater(1.55555, output='SBATCH', input='sec', sround=ceiling)
#' timeFormater(1.55555, output='SBATCH', input='sec', sround=\(x) round(x, 3))
#'
#'
timeFormater <- function(time, output='sec', input = 'min', sround=floor){
    if(output == 'SBATCH'){
        stopifnot(is.numeric(time))
        return(time2SBATCH(time, input=input, sround=sround))
    }
    if(!is.character(time)){
        if(input == 'sec') time <- time / 60
        if(input == 'hour') time <- time * 60
        if(input == 'day') time <- time * 60 * 24
        time <- as.character(time)
    }
    stopifnot(length(time) == 1L && length(output) == 1L)
    stopifnot(output %in% c('sec', 'min', 'hour', 'day'))
    time <- sbatch_time2sec(time)
    if(output == 'min') time <- time / 60
    if(output == 'hour') time <- time / 60 / 60
    if(output == 'day') time <- time / 60 / 60 / 24
    time
}

time2SBATCH <- function(time, input, sround){
    seconds <- switch(input,
                   'day'=time*86400,
                   'hour'=time*3600,
                   'min'=time*60,
                   'sec'=time)
    days <- floor(seconds/86400)
    remainder <- seconds - days*86400
    hours <- floor(remainder / 3600)
    remainder <- remainder - hours * 3600
    minutes <- floor(remainder / 60)
    seconds <- as.character(sround(remainder - minutes * 60))
    minutes <- as.character(minutes)
    if(nchar(seconds) == 1) seconds <- paste0('0', seconds)
    if(nchar(minutes) == 1) minutes <- paste0('0', minutes)
    ret <- sprintf('%s:%s:%s', hours, minutes, seconds)
    if(days > 0)
        ret <- paste0(days, '-', ret)
    ret
}

sbatch_time2sec <- function(time){
    ret <- if(is.character(time)){
        time <- gsub(pattern = " ", "", time)
        time_vec <- c(days=0, hours=0, mins=0, secs=0)
        if(grepl("-", time)){ # day format
            splt <- strsplit(time, "-")[[1L]]
            time_vec['days'] <- as.numeric(splt[1L])
            time <- splt[2L]
            splt <- as.numeric(strsplit(time, ":")[[1L]])
            time <- if(length(splt) == 1L){
                sprintf("%f:00:00", splt[1L])
            } else if(length(splt) == 2L){
                sprintf("%f:%f:00", splt[1L], splt[2L])
            } else if(length(splt) == 3L)
                sprintf("%f:%f:%f", splt[1L], splt[2L], splt[3])
        }
        splt <- as.numeric(strsplit(time, ":")[[1L]])
        time <- if(length(splt) == 1L){
            sprintf("00:%f:00", splt[1L])
        } else if(length(splt) == 2L){
            sprintf("00:%f:%f", splt[1L], splt[2L])
        } else if(length(splt) == 3L){
            time
        } else stop('max_time not correctly specified. Please fix!',
                    call.=FALSE)
        splt <- as.numeric(strsplit(time, ":")[[1L]])
        time_vec[2L:4L] <- splt
        sum(c(86400, 3600, 60, 1) * time_vec)   # c(24*60*60, 60*60, 60, 1)
    } else {
        time
    }
    ret
}

valid_control.list <- function()
    c("stop_on_fatal", "warnings_as_errors", "save_seeds", "store_Random.seeds",
      "store_warning_seeds", "include_replication_index", "include_reps", "try_all_analyse",
      "allow_na", "allow_nan", "type", "print_RAM", "max_time", "max_RAM",
      "tol", "summarise.reg_data", "rel.tol", "k.success", "interpolate.R", "bolster",
      "include_reps", 'global_fun_level', 'useAnalyseHandler', 'max_time.start', 'logging',
      'use_mirai', 'use_load_balancing')

valid_save_details.list <- function()
    c("safe", "compname", "out_rootdir", "save_results_dirname", "save_results_filename",
      "save_seeds_dirname", 'arrayID', "tmpfilename")

on_HPC.cluster <- function() {
    env_vars <- names(Sys.getenv())

    # Common job scheduler prefixes
    scheduler_prefixes <- c(
        "SLURM_",  # Slurm Workload Manager
        "PBS_",    # PBS / Torque
        "LSF_",    # IBM Spectrum LSF (e.g., LSB_JOBID)
        "SGE_"     # Sun Grid Engine / Oracle Grid Engine
    )

    has_scheduler_prefix <- any(sapply(scheduler_prefixes, function(p) any(grepl(paste0("^", p), env_vars))))
    has_lsf_var <- "LSB_JOBID" %in% env_vars
    has_grid_env <- "ENVIRONMENT" %in% env_vars && Sys.getenv("ENVIRONMENT") == "BATCH"

    # Return TRUE if any cluster indicators are found
    return(has_scheduler_prefix || has_lsf_var || has_grid_env)
}

# Test cases:
#
# sbatch_RAM2bytes("1024MB")
# sbatch_RAM2bytes("4G")
# sbatch_RAM2bytes("1.5TB")

sbatch_RAM2bytes <- function(RAM){
    ret <- if(is.character(RAM)){
        RAM <- gsub(pattern = " ", "", RAM)
        type <- logical(3L)
        type[1L] <- grepl('M', RAM)
        type[2L] <- grepl('G', RAM)
        type[3L] <- grepl('T', RAM)
        if(!any(type)) stop('RAM metric must be MB, GB, or TB', call.=FALSE)
        RAM <- gsub('B', "", RAM)
        RAM <- as.numeric(gsub('M|G|T', "", RAM))
        C <- 1000000 # MB2bytes
        if(type[2L]) C <- C * 1000
        if(type[3L]) C <- C * 1000000
        RAM * C
    } else RAM
    ret
}

#' @rdname genSeeds
#' @param ... does nothing
gen_seeds <- function(...){
    .Deprecated('genSeeds')
    genSeeds(...)

}

add_cbind <- function(lst){
    len <- sapply(lst, ncol)
    if(!any(len)) return(lst[[1L]])
    lst <- lapply(lst, \(x){
        x[is.na(x)] <- 0
        x
    })
    for(i in 1L:length(lst)){
        if(length(lst[[i]])){
            ret <- lst[[i]]
            if(i == length(lst)) return(ret)
            break
        }
    }
    from <- i + 1L
    for(i in from:length(lst)){
        nms <- colnames(ret)
        nms2 <- colnames(lst[[i]])
        matched <- nms %in% nms2
        if(any(matched)){
            for(j in 1L:length(nms))
                if(matched[j])
                    ret[,nms[j]] <- ret[,nms[j]] + lst[[i]][,nms[j]]
        }
        ret <- cbind(ret, lst[[i]][,!(nms2 %in% nms)])
    }
    dplyr::as_tibble(ret)
}


#' Read simulation files
#'
#' Convenience function that switches between \code{readRDS()} and
#' \code{qs2::qs_read()} for functions saved with \code{SimDesign}.
#' By convention, objects saved with the extension \code{.rds} are read
#' as R binary files and typically reflect the final object from
#' \code{\link{runSimulation}} or \code{\link{runArraySimulation}},
#' while files without an extension are read using \code{qs2} (most often
#' temporary files or results written to an associated sub-directory).
#'
#' @export
#' @return the R binary object
#' @param filename name of the file to read in
SimRead <- function(filename){
    file_ext <- tools::file_ext(filename)
    tmp <- if(tolower(file_ext) == 'rds')
        try(readRDS(filename), TRUE)
    else try(qs2::qd_read(filename), TRUE)
    tmp
}
