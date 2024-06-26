% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{Generate}
\alias{Generate}
\title{Generate data}
\usage{
Generate(condition, fixed_objects)
}
\arguments{
\item{condition}{a single row from the \code{design} input (as a \code{data.frame}), indicating the
simulation conditions}

\item{fixed_objects}{object passed down from \code{\link{runSimulation}}}
}
\value{
returns a single object containing the data to be analyzed (usually a
  \code{vector}, \code{matrix}, or \code{data.frame}),
  or \code{list}
}
\description{
Generate data from a single row in the \code{design} input (see \code{\link{runSimulation}}). R contains
numerous approaches to generate data, some of which are contained in the base package, as well
as in \code{SimDesign} (e.g., \code{\link{rmgh}}, \code{\link{rValeMaurelli}}, \code{\link{rHeadrick}}).
However the majority can be found in external packages. See CRAN's list of possible distributions here:
\url{https://CRAN.R-project.org/view=Distributions}. Note that this function technically
can be omitted if the data generation is provided in the \code{\link{Analyse}} step, though
in general this is not recommended.
}
\details{
The use of \code{\link{try}} functions is generally not required in this function because \code{Generate}
is internally wrapped in a \code{\link{try}} call. Therefore, if a function stops early
then this will cause the function to halt internally, the message which triggered the \code{\link{stop}}
will be recorded, and \code{Generate} will be called again to obtain a different dataset.
That said, it may be useful for users to throw their own \code{\link{stop}} commands if the data
should be re-drawn for other reasons (e.g., an estimated model terminated correctly
but the maximum number of iterations were reached).
}
\examples{
\dontrun{

generate <- function(condition, fixed_objects) {
    N1 <- condition$sample_sizes_group1
    N2 <- condition$sample_sizes_group2
    sd <- condition$standard_deviations

    group1 <- rnorm(N1)
    group2 <- rnorm(N2, sd=sd)
    dat <- data.frame(group = c(rep('g1', N1), rep('g2', N2)),
                      DV = c(group1, group2))
    # just a silly example of a simulated parameter
    pars <- list(random_number = rnorm(1))

    list(dat=dat, parameters=pars)
}

# similar to above, but using the Attach() function instead of indexing
generate <- function(condition, fixed_objects) {
    Attach(condition)
    N1 <- sample_sizes_group1
    N2 <- sample_sizes_group2
    sd <- standard_deviations

    group1 <- rnorm(N1)
    group2 <- rnorm(N2, sd=sd)
    dat <- data.frame(group = c(rep('g1', N1), rep('g2', N2)),
                      DV = c(group1, group2))
    dat
}

generate2 <- function(condition, fixed_objects) {
    mu <- sample(c(-1,0,1), 1)
    dat <- rnorm(100, mu)
    dat        #return simple vector (discard mu information)
}

generate3 <- function(condition, fixed_objects) {
    mu <- sample(c(-1,0,1), 1)
    dat <- data.frame(DV = rnorm(100, mu))
    dat
}

}

}
\references{
Chalmers, R. P., & Adkins, M. C.  (2020). Writing Effective and Reliable Monte Carlo Simulations
with the SimDesign Package. \code{The Quantitative Methods for Psychology, 16}(4), 248-280.
\doi{10.20982/tqmp.16.4.p248}

Sigal, M. J., & Chalmers, R. P. (2016). Play it again: Teaching statistics with Monte
Carlo simulation. \code{Journal of Statistics Education, 24}(3), 136-156.
\doi{10.1080/10691898.2016.1246953}
}
\seealso{
\code{\link{addMissing}}, \code{\link{Attach}},
  \code{\link{rmgh}}, \code{\link{rValeMaurelli}}, \code{\link{rHeadrick}}
}
