\name{pheno}
\alias{pheno}

\title{
  Generates posterior samples from a phenology model
}

\description{
  Generates posterior samples from a phenology model using an option of
  several different likelihoods and prior configurations.
}

\usage{
pheno(formula, data = parent.frame(), gamma = c(0, 1), family = "normal",
      starting, tuning, priors, n.samples, sub.sample,
      summary.only = FALSE, fitted = FALSE, verbose = TRUE, n.report=100, ...)
}

\arguments{
  \item{formula}{a symbolic description of the regression model to be
  fit. The left hand size is the name of the Vegetation Index (VI) variable and the
  right and size is the name of the time variable.}

  \item{data}{an optional data frame containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{pheno} is called. }

  \item{gamma}{a vector of length two that holds the VI
  variable's theoretical bounds. Lower and upper bounds are given in
  element 1 and 2, respectively. These must be finite values.}
  
  \item{family}{a quoted string that indicates which likelihood to use
  for modeling the VI variable. Options are Beta, Normal, and
  truncated Normal, specified using quoted argument values
  \code{beta}, \code{normal}, and \code{t.normal}, respectively. If
  \code{family} equals \code{"beta"} then The \eqn{\sigma^2} parameter
  is the variance for the Beta, Normal, and truncated Normal likelihoods. See Details for VI variable's support for each likelihood.}

  \item{starting}{a list with each tag corresponding to a parameter
  name. Valid tags are \code{alpha.1}, \code{alpha.2}, \ldots,
  \code{alpha.7}, and \code{sigma.sq}.  The value portion of each tag is the
  parameter's starting value. If \code{starting} is not provided, then each parameter's
  starting value is set as a random value from its prior distribution (in
  most cases this does not provide good starting values). If
  \eqn{\alpha} starting values are provided, they must be within their
  respective prior distributions, see the Details for additional
  guidance on selecting starting values. }

  \item{tuning}{a list with each tag corresponding to a parameter
  name. Valid tags are \code{alpha.1}, \code{alpha.2}, \ldots,
  \code{alpha.7}, and \code{sigma.sq}.  The value portion of each tag defines
  the variance of the Metropolis sampler Normal proposal
  distribution. Tuning values should be selected to keep the acceptance
  rate between approximately 20 and 50 percent.}

  \item{priors}{a list with tags \code{alpha} and
  \code{sigma.sq.IG}. \code{alpha} should be a list comprising
  \code{alpha.1.Unif}, \code{alpha.2.Unif}, \ldots,
  \code{alpha.7.Unif} with each of these tags set equal to a vector of
  length two that holds the lower and upper bound for the Uniform  priors on the given \code{alpha}. The \code{sigma.sq.IG} is a vector of
  length two that holds the shape and scale parameters for the
  \eqn{\sigma^2}'s inverse-Gamma prior distribution. At minimum the \code{priors} list
  must include a prior for \code{sigma.sq.IG}. If \eqn{\alpha} priors are not specified then they are assumed to follow their default values, see the Details section below.}

  \item{n.samples}{the number of MCMC samples to collect.}

  \item{sub.sample}{an optional list that specifies the subset of MCMC samples should be returned. Valid tags are \code{start}, \code{end}, and \code{thin}. The default values are \code{start=floor(0.5*n.samples)}, \code{end=n.samples} and \code{thin=1}.}

  \item{fitted}{if \code{TRUE}, regression fitted are returned. The
  argument \code{sub.sample} controls which MCMC samples are used to
  generate the fitted values.}

  \item{summary.only}{if \code{TRUE}, returns the mean, standard deviation, and quantiles of MCMC samples as well as fitted values if \code{fitted=TRUE}. The summaries are based on the entire chain or subsamples if \code{sub.sample} is provided. The defualt probabilities
  used in the \code{quantile} function are 0.025, 0.5, and 0.975. These
  default quantiles can be changed, see Details.}

  \item{verbose}{if \code{TRUE}, model specification and progress of the
  sampler is printed to the screen. Otherwise, nothing is printed to the screen.}

  \item{n.report}{the interval to report Metropolis sampler acceptance and MCMC progress.}

  \item{...}{optional argument \code{summary.prob} is used to pass a
    vector of probabilities to the \code{quantile} function's
    \code{probs} argument which is used to compute summaries of MCMC samples
    and fitted values if \code{fitted = TRUE}. Optional argument
    \code{t.normal.bounds} takes a vector of length two with the
    first and second elements defining the truncated Normal's lower and upper bounds, respectively. The default bounds for the truncated
    Normal are (0,1).} 
}

\details{
  
Selection of the likelihood via the \code{family} argument should be
respective of possible VI value range. The Beta likelihood
assumes support for VI between 0 and 1. The Normal likelihood assumes support for VI on
the whole real line. The default for the truncated Normal likelihood is
support between 0 and 1; however, specifying the optional argument
\code{t.normal.bounds} allows for user defined support bounds.
  
The default priors for the \eqn{\alpha}'s are:
\tabular{lll}{
  \eqn{\alpha_1 \sim Unif(\gamma_1, \gamma_2)}{alpha.1 ~ Unif(gamma[1], gamma[2])},\tab\eqn{\alpha_2 \sim Unif(0, \gamma_2-\alpha_1)}{alpha.2 ~ Unif(0, gamma[2]-alpha.1)},\tab\eqn{\alpha_3 \sim Unif(0, 1)}{alpha.3 ~ Unif(0, 1)}\cr
  \eqn{\alpha_4 \sim Unif(0, \alpha_7)}{alpha.4 ~ Unif(0, alpha.7)},\tab\eqn{\alpha_5 \sim Unif(0, 0.001)}{alpha.5 ~ Unif(0, 0.001)},\tab\eqn{\alpha_6 \sim Unif(0, 1)}{alpha.6 ~ Unif(0, 1)}\cr
  \eqn{\alpha_7 \sim Unif(1, 365)}{alpha.7 ~ Unif(0, 1)}\tab\tab\cr
}
The hyperparameters of these prior distributions can be change using the
\code{priors} argument.

An error will be thrown if an \eqn{\alpha}{alpha} starting value
provided in the \code{starting} list is outside the default or supplied
priors. Note, \eqn{\alpha_2}{alpha.2} and \eqn{\alpha_4}{alpha.4} have
upper bounds determined by starting values of \eqn{\alpha_1}{alpha.1}
and \eqn{\alpha_7}{alpha.7}, respectively. As noted in the description
of the \code{starting} argument, starting values will be automatically
generated if not provided. We strongly suggest that reasonable starting
values are provided or the MCMC sampler might not perform well.
}
\value{
The output is a list of class \code{pheno}. List components include:
\item{p.theta.samples}{is a \code{coda} object of \code{n.samples} (or
  number of subsamples specified in the optional argument
  \code{sub.samples}) with columns labeled \code{alpha.1}, \code{alpha.2}, \ldots,
  \code{alpha.7}, and \code{sigma^2}.}
\item{MH.acceptance}{is a named vector with values giving the Metropolis
  samplers \code{overall} acceptance percent and the \code{last batch}
  which is the percent acceptance in the last \code{n.report} samples.}
\item{p.fitted.samples}{if \code{fitted=TRUE} then this \eqn{n}, i.e.,
  the number of VI observations, by \code{n.samples} (or
  number of subsamples specified in the optional argument
  \code{sub.samples}) matrix holds the posterior fitted values sampled
  one-for-one from the parameters posterior values.} 
}

\examples{
  
##select a single pixel and a few years
df <- subset(aoi_hls, pixel == 1 & year \%in\% c(2017,2018,2019))
plot(df$doy, df$evi)

set.seed(1)

m.beta <- pheno(evi~doy, data=df, family="beta",
                starting=list(alpha.1=0.1, alpha.2=0.5, alpha.3=0.25,
                         alpha.4=50, alpha.5=0.0001, alpha.6 = 0.25, alpha.7=200, sigma.sq=0.01),
		tuning=list(alpha.1=0.001, alpha.2=0.01, alpha.3=0.001,
			 alpha.4=0.1, alpha.5=0.00005, alpha.6=0.001, alpha.7=0.5, sigma.sq=0.1),
                n.samples=50000, n.report=25000,
                priors=list(alpha=list(alpha.5=c(-1, 1)), sigma.sq.IG=c(2, 0.01)),
                fitted = TRUE, sub.sample=list("start"=25000)
                )

summary(m.beta)

}

\author{
  Andrew O. Finley \email{finleya@msu.edu},
  Chad Babcock \email{cbabcock@umn.edu} 
}
