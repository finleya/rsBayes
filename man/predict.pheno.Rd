\name{predict.pheno}

\alias{predict.pheno}

\title{Function for prediction at new time points using \code{pheno} models.}

\description{
  The function \code{predict} collects posterior predictive samples
  for a set of new time points given an object of
  class \code{pheno}.
}

\usage{
\method{predict}{pheno}(object, t.0, sub.sample, ...)
}

\arguments{
  \item{object}{an object of class \code{pheno}. }
  \item{t.0}{the vector of time points for which posterior predictive
    samples will be generated. }
  \item{sub.sample}{an optional list that specifies the samples to included in
    the composition sampling a non-Conjugate model. Valid tags are \code{start},
    \code{end}, and \code{thin}. Given the value associated with the tags,
    the sample subset is selected using \code{seq(as.integer(start),
      as.integer(end), by=as.integer(thin))}. The default values are
    \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and
    \code{thin=1}.}  
  \item{...}{See details for additional optional arguments.}
}

\details{
  The following optional arguments can be passed:
  
  \code{n.omp.threads} which is a positive integer indicating
  the number of threads to use for SMP parallel processing. The package must
  be compiled for OpenMP support. For most Intel-based machines, we recommend setting
  \code{n.omp.threads} up to the number of hyperthreaded cores. Note, \code{n.omp.threads} > 1 might not
  work on some systems.
  
  \code{verbose} which if \code{TRUE}, the progress of the
  sampler is printed to the screen. Otherwise, nothing is printed to
  the screen.  
  
  \code{n.report} which is the interval to report sampling progress. 
}

\value{
  An object of class \code{predict.pheno} which is a list comprising:
  
  \item{p.predictive.samples}{a matrix that holds the response variable posterior
    predictive samples where rows are time points corresponding to \code{t.0}.}
  \item{p.predictive.quantiles}{a matrix that holds quantiles of the response variable posterior
    predictive samples where rows are time points corresponding to \code{t.0}.}
  \item{run.time}{execution time reported using \code{proc.time()}.}
}

\author{
  Andrew O. Finley \email{finleya@msu.edu}, \cr
  Chad Babcock \email{cbabcock@umn.edu}
}

\keyword{model}