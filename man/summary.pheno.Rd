\name{summary.pheno}

\alias{summary.pheno}
\alias{plot.pheno}
\alias{fitted.pheno}
\alias{residuals.pheno}
\alias{print.pheno}
\alias{print.predict.pheno}

\title{Methods for pheno and derived objects}

\description{
  Methods for extracting information from fitted model of class
  \code{pheno} and \code{predict.pheno} objects from \code{predict}.
}

\usage{
\method{summary}{pheno}(object, sub.sample, quantiles = c(0.025, 0.25,
0.5, 0.75, 0.975), digits = max(3L, getOption("digits") - 3L), \dots)
\method{plot}{pheno}(x, \dots)
\method{fitted}{pheno}(object, sub.sample, \dots)
\method{residuals}{pheno}(object, sub.sample, \dots)
\method{print}{pheno}(x, \dots)
\method{print}{predict.pheno}(x, \dots)
}

\arguments{
  \item{object, x}{object of class \code{pheno} or \code{predict.pheno}.}
  \item{sub.sample}{an optional list that specifies the samples to included in
    the summary or composition sampling. Valid tags are \code{start},
    \code{end}, and \code{thin}. Given the value associated with the tags,
    the sample subset is selected using \code{seq(as.integer(start),
      as.integer(end), by=as.integer(thin))}. The default values are
    \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and \code{thin=1}.}
  \item{quantiles}{for \code{summary}, posterior distribution quantiles to compute.}
  \item{digits}{for \code{summary}, number of digits to report in summary.}
  \item{\dots}{currently not used.}
}

\details{
  A set of standard extractor functions for fitted model objects of
  class \code{pheno} and prediction object of class \code{predict.pheno}, including methods to the generic functions \code{\link[base]{print}} and \code{\link[base]{summary}}.
}

\author{
  Andrew O. Finley \email{finleya@msu.edu},
  Chad Babcock \email{cbabcock@umn.edu} 
}

\keyword{model}