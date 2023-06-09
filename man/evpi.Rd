% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evpi.R
\name{evpi}
\alias{evpi}
\title{Calculate the expected value of perfect information from a decision model}
\usage{
evpi(outputs, nsim = NULL)
}
\arguments{
\item{outputs}{This could take one of two forms

"net benefit" form: a matrix or data frame of samples from the uncertainty
distribution of the expected net benefit.  The number of rows should equal
the number of samples, and the number of columns should equal the number
of decision options.

"cost-effectiveness analysis" form: a list with the following named
components:

\code{"c"}: a matrix or data frame of samples from the distribution of
costs.  There should be one column for each decision option.

\code{"e"}: a matrix or data frame of samples from the distribution of
effects, likewise.

\code{"k"}: a vector of willingness-to-pay values.

Objects of class \code{"bcea"}, as created by the \pkg{BCEA} package, are in
this "cost-effectiveness analysis" format, therefore they may be supplied as
the \code{outputs} argument.

Users of \pkg{heemod} can create an object of this form, given an object
produced by \code{run_psa} (\code{obj}, say), with \code{\link{import_heemod_outputs}}.

If \code{outputs} is a matrix or data frame, it is assumed to be of "net
benefit" form.  Otherwise if it is a list, it is assumed to be of "cost
effectiveness analysis" form.}

\item{nsim}{Number of simulations from the decision model to use
for calculating EVPPI.  The first \code{nsim} rows of the
objects in \code{inputs} and \code{outputs} are used.}
}
\value{
The expected value of perfect information, either as a single value, or a data frame indicating the value for each willingness-to-pay.
}
\description{
Calculate the expected value of perfect information from a decision model using standard Monte Carlo simulation
}
