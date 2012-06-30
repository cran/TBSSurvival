% File TBSSurvival/man/tbs.survreg.be.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{tbs.survreg.be}

%%%
\alias{tbs.survreg.be}

\title{Bayesian Estimation of the TBS Model for Failure Data}
\description{
  This function perform the Bayesian estimation of the TBS model. The prior for 
  the parameters `lambda' and `xi' are uniform-exponential mixture and, if not 
  specified, for parameter beta0 is a normal with mean 5 and sd 5. The estimations are
  done by metroplois-hasting (using the function `metrop' availible at package `mcmc').
}
\usage{
  tbs.survreg.be(formula, dist = "norm",  max.time = -1, guess.beta, guess.lambda,
                 guess.xi, burn = 1000, jump = 2, size = 1000, scale = 1, 
                 prior.mean = NULL, prior.sd = NULL)
}
\arguments{
  \item{formula}{A formula specification containing a Surv model with right-censored (or no censored) data as in the package survival.}
  \item{dist}{error distribution; \code{dist} = \code{"norm"}, \code{"doubexp"}, \code{"t"}, \code{"cauchy"} or \code{"logistic"}.}
  \item{max.time}{Maximum time (in minutes) to run the optimization (<= 0 means no limit).}
  \item{guess.beta}{initial value of the Markov Chain for the vector `beta'.}
  \item{guess.lambda}{initial value of the Markov Chain for the parameter `lambda'.}
  \item{guess.xi}{initial value of the Markov Chain for the parameter `xi'.}
  \item{burn}{burn-in, number of firsts samples of posterior to not use.}
  \item{jump}{number of jump between each sample of posterior to avoid the problem of auto-correlation
              between the samples.}
  \item{size}{size of final sample of posterior.}
  \item{scale}{parameter of `metrop' function. Controls the acceptance rate.}
  \item{prior.mean}{Prior Mean for the MCMC.}
  \item{prior.sd}{Prior std deviation for the MCMC.}
}
\details{
  This function was developed to the problem without co-variables. Some changes maybe are need to work
  with co-variables.
}
\value{
  A list with components:
  \item{call}{function evaluated.}
  \item{x}{co-variable matrix used.}
  \item{time}{survival time.}
  \item{delta}{censor status.}
  \item{post}{posterior sample of the parameters.}
  \item{par}{the parameters are estimated as median of posterior sample.}
  \item{par.sd}{standard deviation of posterior sample.}
  \item{par.HPD}{95\% high posterior density credal interval of each parameter.}
  \item{DIC}{Deviance Information Criterion.}
  \item{error}{summary statistics and high posterior density credal interval of 95\% for
                  the posterior of error of TBS model.}
  \item{run.time}{Time spent with the function.}
}
\references{
  Meeker, W. and Escobar, L. (1998) \emph{Statistical Methods for Reliability Data}. Willey, ISBN 0-471-14328-6.
}
\seealso{
  \code{\link{tbs.survreg.mle},\link{dtbs},\link{ptbs},\link{qtbs},\link{rtbs}}.
}
\examples{
# set.seed is used to produce the same results all times.
set.seed(1234)

# Alloy - T7987: data extracted from Meeker and Escobar (1998), pp. 131)
data(alloyT7987)
alloyT7987$time  <- as.double(alloyT7987$time)
alloyT7987$delta <- as.double(alloyT7987$delta)

# Bayesian estimation with logistic error
formula <- Surv(alloyT7987$time,alloyT7987$delta == 1) ~ 1
tbs.be <- tbs.survreg.be(formula,guess.lambda=1,guess.xi=1,guess.beta=5,
                         dist="logistic",burn=1000,jump=10,size=1000,scale=0.06)
# evaluating the estimated survival and HPD
axis.x <- seq(0.01,300,1)
aux <- matrix(NA,length(axis.x),length(tbs.be$post[,1]))
for (j in 1:length(tbs.be$post[,1])) {
  aux[,j] <- 1-ptbs(axis.x,lambda=tbs.be$post[j,1],xi=tbs.be$post[j,2],
                    beta=tbs.be$post[j,3],dist="logistic")
}
survival <- matrix(NA,length(axis.x),3)
for (i in 1:length(axis.x)) {
  survival[i,] <- c(mean(aux[i,]),HPDinterval(as.mcmc(aux[i,]),0.95))
}
rm(aux,i,j)

# Kapan-Meier estimator
km <- survfit(formula = Surv(alloyT7987$time, alloyT7987$delta == 1) ~ 1)

plot(km,ylab="",xlab="", xlim=c(min(alloyT7987$time),max(alloyT7987$time)),
     conf.int=FALSE,axes=FALSE,lty=1,lwd=1)
t <- seq(min(alloyT7987$time),max(alloyT7987$time),
         (max(alloyT7987$time)-min(alloyT7987$time)-0.01)/1000)
title(ylab="R(t)", xlab="t: number of cycles (in thousands)",
      sub="Alloy - T7987 (cf. Meeker and Escobar (1998), pp. 131)",
      main="Reliability function (BE)",cex.lab=1.2)
axis(1,at=c(93,100,150,200,250,300),lwd=2,lwd.ticks=2,pos=0)
axis(2,lwd=2,lwd.ticks=2,pos=min(alloyT7987$time))
legend(170,0.95,c("Kaplan-Meier",
                  expression(textstyle(paste("TBS / ",sep="")) ~ epsilon
                  ~ textstyle(paste("~",sep="")) ~ Logistic),
                  "0.95 HPD Interval"),
       col=c(1,2,2),lty=c(1,1,2),cex=1.1,lwd=c(1,2,2),bg="white")
lines(axis.x,survival[,1],type="l",lwd=2,col=2,lty=1)
lines(axis.x,survival[,2],type="l",lwd=2,col=2,lty=2)
lines(axis.x,survival[,3],type="l",lwd=2,col=2,lty=2)
}
