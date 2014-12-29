% File TBSSurvival/man/dweib.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012-2013, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                          Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{tbs.survreg.mle}

%%%
\alias{tbs.survreg.mle}

\title{MLE of the TBS Model for Failure Data}
\description{
  This function performs the Maximum Likelihhod Estimation of the TBS model.
  The optimization is done by the function `optim' (or optionally the package Rsolnp when available).
}
\usage{
 tbs.survreg.mle(formula,dist=dist.error("all"),
                 method=c("Nelder-Mead", "BFGS", "Rsolnp", "SANN", "CG"),
                 verbose=FALSE,nstart=10,max.time=-1,seed=1234,gradient=FALSE)
}
\arguments{
  \item{formula}{A formula specification containing a Surv model with right-censored data as in the package survival.}
  \item{dist}{error distribution; \code{dist} can be given by name (\code{"norm"}, \code{"doubexp"}, \code{"t"}, \code{"cauchy"} or \code{"logistic"}) or by \link{dist.error}.}
  \item{method}{a vector of numerical methods to be used in the optimization. The function try all listed methods and returns all results, together with an indication of the
                solution with maximal likelihood among them.}
  \item{verbose}{Boolean to indicate the amount of output during the execution of the optimization.}
  \item{nstart}{Number of feasible initial points to guess when performing the optimization.}
  \item{max.time}{Maximum time (in minutes) to run the optimization (<= 0 means no limit).}
  \item{seed}{The number that is used to initialize the seed for random number generation.}
  \item{gradient}{if TRUE, MLE tries to use the implemented gradient functions (usually the numerical ones are ok).}
}
\details{
  This function calls numerical optimization methods to maximize the likelihood 
  of the TBS model, according to the given error distribution, method of optimization, 
  and formula. The formula is supposed to have a Surv object and possibility co-variates,
  just as the standard specification of R formulas. The optimizers are going to do their 
  best to find high likelihood estimates, but as in most estimation methods that need 
  a numerical optimization procedure, the obtained estimate cannot be guaranteed to be a
  global optimal solution, but instead is dependent on the initial guessing points, and thus on 
  the seed of the random number generation.
}
\value{
  Either an element of class tbs.survreg.mle (with print, summary, and plot functions) or a list of them (depending whether the
  call of tbs.survreg.mle was made for a single distribution or a list of them). In case it is a list, additional fields named
  best and best.n give the name and the position of the best estimation in the list, respectively. Each element of tbs.survreg.mle has the following components:
  \item{lambda}{The estimate for parameter lambda}
  \item{xi}{The estimate for parameter xi}
  \item{beta}{A vector with the estimate for parameter beta}
  \item{lambda.se}{The standard error for parameter lambda}
  \item{xi.se}{The standard error for parameter xi}
  \item{beta.se}{A vector with the standard error for parameter beta}
  \item{log.lik}{The log-likelihood at parameters \code{par}.}
  \item{error.dist}{The error distribution chosen.}
  \item{AIC}{Akaike Information Criterion.}
  \item{AICc}{AICc is AIC with a second order correction for small sample sizes.}
  \item{BIC}{Bayesian Information Criterion.}
  \item{method}{Numerical method used to achive the MLE.}
  \item{convergence}{If \code{convergence} is \code{FALSE} then it was not possible to find the MLE.}
  \item{time}{observed survival times.}
  \item{error}{error of the estimated model.}
  \item{call}{function evaluated.}
  \item{formula}{formula entered by user.}
  \item{run.time}{Time spent with the function.}
}
\references{
  Meeker, W. and Escobar, L. (1998) \emph{Statistical Methods for Reliability Data}. Willey, ISBN 0-471-14328-6.
}
\seealso{
  \code{\link{dist.error},\link{tbs.survreg.be},\link{dtbs},\link{ptbs},\link{qtbs},\link{rtbs}}.
}
\examples{
# Alloy - T7987: data extracted from Meeker and Escobar (1998), pp. 131.
data(alloyT7987)
alloyT7987$time <- as.double(alloyT7987$time)
alloyT7987$delta <- as.double(alloyT7987$delta)

# MLE estimation with logistic error
formula <- survival::Surv(alloyT7987$time,alloyT7987$delta == 1) ~ 1
tbs.mle <- tbs.survreg.mle(formula,dist=dist.error("logistic"),method="Nelder-Mead",nstart=3)

# Kaplan-Meier estimation
km <- survival::survfit(formula)

# Plot survival function
plot(tbs.mle)
lines(km)

# Plot hazard function
plot(tbs.mle, plot.type="hazard")

}

