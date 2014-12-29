% File TBSSurvival/man/tbs.survreg.be.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012-2013, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{tbs.survreg.be}

%%%
\alias{tbs.survreg.be}

\title{Bayesian Estimation of the TBS Model for Survival Data}
\description{
  This function performs the Bayesian estimation of the Transform-Both-Sides (TBS) model. The priors for 
  the parameters `lambda' and `xi' are uniform-exponential mixtures and, if not 
  specified, for parameter beta is a normal with mean 5 and sd 5. The estimations are
  done by Metropolis-Hasting (using the function `metrop' availible with the package `mcmc').
}
\usage{
  tbs.survreg.be(formula, dist=dist.error("norm"),max.time = -1, guess.beta = NULL, 
                 guess.lambda = 1, guess.xi = 1, burn = 1000, jump = 2, size = 500, 
                 scale = 0.1, prior.mean = NULL, prior.sd = NULL, seed = 1234)
}
\arguments{
  \item{formula}{A formula specification containing a \link{Surv} model with right-censored (or no censored) data as in the package survival.}
  \item{dist}{Error distribution; \code{dist} can be given by name (\code{"norm"}, \code{"doubexp"}, \code{"t"}, \code{"cauchy"} or \code{"logistic"}) or by \link{dist.error}.}
  \item{max.time}{Maximum time (in minutes) to run the optimization (<= 0 means no limit).}
  \item{guess.beta}{Initial value of the Markov Chain for the vector `beta'. Default will fill it with zeros.}
  \item{guess.lambda}{Initial value of the Markov Chain for the parameter `lambda'.}
  \item{guess.xi}{Initial value of the Markov Chain for the parameter `xi'.}
  \item{burn}{Burn-in: number of initial samples of the posterior not to use.}
  \item{jump}{Number of jumps between each sample of the posterior to avoid the problem of auto-correlation
              between the samples.}
  \item{size}{Size of final sample of the posterior.}
  \item{scale}{Parameter of `metrop' function. Controls the acceptance rate.}
  \item{prior.mean}{Prior Mean for the MCMC.}
  \item{prior.sd}{Prior std deviation for the MCMC.}
  \item{seed}{The number that is used to initialize the seed for random number generation.}
}
\details{
  This function performs the Bayesian estimation of the Transform-Both-Sides (TBS) model. The priors for 
  the parameters `lambda' and `xi' are uniform-exponential mixtures and, if not 
  specified, for parameter beta is a normal with mean 5 and sd 5. The estimations are
  done by Metropolis-Hasting (using the function `metrop' availible with the package `mcmc').
}
\value{
  An element of the class tbs.survreg.be, with the components:
  \item{call}{function evaluated.}
  \item{x}{co-variable matrix used.}
  \item{time}{survival time.}
  \item{delta}{censor status.}
  \item{post}{posterior sample of the parameters.}
  \item{lambda}{posterior mean of lambda.}
  \item{xi}{posterior mean of xi.}
  \item{beta}{vector with posterior mean of beta.}
  \item{lamda.sd}{standard deviation for lambda.}
  \item{xi.sd}{standard deviation of for xi.}
  \item{beta.sd}{standard deviation of for beta.}
  \item{lambda.HPD}{95\% high posterior density credal interval of lambda.}
  \item{xi.HPD}{95\% high posterior density credal interval of xi.}
  \item{beta.HPD}{95\% high posterior density credal interval vector of beta.}
  \item{DIC}{Deviance Information Criterion.}
  \item{error}{summary statistics for the posterior of error of TBS model.}
  \item{error.dist}{error distribution.}
  \item{run.time}{Time spent with the function.}
}
\references{
  Meeker, W. and Escobar, L. (1998) \emph{Statistical Methods for Reliability Data}. Willey, ISBN 0-471-14328-6.
}
\seealso{
  \code{\link{dist.error},\link{tbs.survreg.mle},\link{dtbs},\link{ptbs},\link{qtbs},\link{rtbs}}.
}
\examples{
# set.seed is used to produce the same results all times.
set.seed(1234)

# Alloy - T7987: data extracted from Meeker and Escobar (1998), pp. 131)
data(alloyT7987)
alloyT7987$time  <- as.double(alloyT7987$time)
alloyT7987$delta <- as.double(alloyT7987$delta)

# Bayesian estimation with logistic error
formula <- survival::Surv(alloyT7987$time,alloyT7987$delta == 1) ~ 1
tbs.be <- tbs.survreg.be(formula,guess.lambda=1,guess.xi=1,guess.beta=5,
                         dist=dist.error("logistic"),burn=1000,jump=10,size=500,scale=0.06)

# Kapan-Meier estimator
km <- survival::survfit(formula = survival::Surv(alloyT7987$time, alloyT7987$delta == 1) ~ 1)

# Plot survival function
plot(tbs.be,lwd=2,HPD=TRUE,HPD.alpha=0.95,col.HPD=2,lty.HPD=1,lwd.HPD=2)
lines(km)

# Plot survival function
plot(tbs.be,plot.type="hazard",lwd=2,HPD=TRUE,HPD.alpha=0.95,col.HPD=2,lty.HPD=1,lwd.HPD=2)

# Plot auto-correlation of the posterior sample
plot(tbs.be,plot.type="auto")

# Plot "time-series" of the posterior sample
plot(tbs.be,plot.type="ts")
}
