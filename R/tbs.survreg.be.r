# TBSSurvival package for R (http://www.R-project.org)
# Copyright (C) 2012 Adriano Polpo, Cassio de Campos, Debajyoti Sinha
#                    Jianchang Lin and Stuart Lipsitz.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

## TBS estimation using a Bayesian approach. The lack of a closed form
## solution forces us to tackle the problem by MCMC.
## max.time is a time limit in minutes (<= 0 means no limit)
## formula is specification containing a Surv model with right-censored data as in the package survival.
## dist defines the error distribution: "norm", "doubexp", "t", "cauchy", "logistic".
## guess.beta, guess.lambda, guess.xi are initial value of the Markov Chain (beta has to have the same
## of elements as covariates, lambda and xi are scalars.
## burn-in is the number of firsts samples of posterior to not use, and jump is the number of jump
## between each sample of posterior to avoid the problem of auto-correlation between the samples.
## size is the final sample size of the posterior. scale is a parameter of the `metrop' function which
## controls the acceptance rate. prior.mean and prior.sd define the parameters of the normal prior for
## the MCMC (by default, they are equal to 5 and 5).
tbs.survreg.be <- function(formula,dist="norm",max.time=-1,
                           guess.beta,guess.lambda,guess.xi,
                           burn=1000,jump=2,size=1000,scale=1,
                           prior.mean=NULL,prior.sd=NULL) {
  initial.time <- .gettime()
  if(max.time <= 0) {
    ## user didn't define a timeout, so we set to a large number
    max.time <- 1e10
  }

  ## check the class of formula
  if (attributes(formula)$class != "formula")
    stop("A formula argument is required")

  ## record the call arguments
  Call  <- match.call()
  ## read the information from within the formula to populate the required variables
  mf <- model.frame(formula=formula)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  time <- y[,1]
  delta <- y[,2]
  x.k   <- dim(x)[2]
  n     <- dim(x)[1]
  ## check if delta is an indicator function
  if (any((delta != 0) & (delta != 1)))  {
    stop("It is only accepted uncensored or right censored data")
  }

  out <- NULL
  out$call <- Call
  if (is.matrix(x))
    out$x <- x[order(time),]
  else
    out$x <- x[order(time)]
  out$delta <- delta[order(time)]
  out$time <- time[order(time)]

  ## perform a series of verifications for the given arguments of the function
  if (length(guess.lambda) != 1)
    stop("guess.lambda is not a scalar")
  if (guess.lambda <= 0)
    stop("guess.lambda must be a positive number")
  if (length(guess.xi) != 1)
    stop("guess.xi is not a scalar")
  if (guess.xi <= 0)
    stop("guess.xi must be a positive number")
  if (length(guess.beta) != x.k)
    stop("guess.beta length is not in accordance with the model specification")
  guess <- c(guess.lambda,guess.xi,guess.beta)

  if ((!is.integer(burn)) && (burn <= 0)) 
    stop("burn must be a positive integer number")
  if ((!is.integer(jump)) && (jump <= 0)) 
    stop("jump must be a positive integer number")
  if ((!is.integer(size)) && (size <= 0)) 
    stop("size must be a positive integer number")
  if (scale <= 0)
    stop("scale must be a positive number")

  if (!is.null(x)) {
    if (is.matrix(x)) {
      if (length(time) != length(x[,1]))
        stop("length of time is different of length of x")
    } else {
      if (length(time) != length(x))
        stop("length of time is different of length of x")
      x <- matrix(x,length(x),1)
    }
    if(length(beta) > 1)
      beta <- matrix(beta,length(beta),1)
  } else {
    x <- matrix(1,length(time),1)
  }

  if (is.null(prior.mean)) {
    prior.mean <- 5
  } else {
    if (!is.vector(prior.mean))
      stop("mean is not a vector/scalar")
    if ((length(prior.mean) != 1) || (length(prior.mean) != length(guess[3:length(guess)]))) {
      stop(paste("length mean is different of 1 or ",length(guess[3:length(guess)]),sep=""))
    }
  }
  if (is.null(prior.sd)) {
    prior.sd <- 5
  } else {
    if (!is.vector(prior.sd))
      stop("sd is not a vector/scalar")
    if ((length(prior.sd) != 1) || (length(prior.sd) != length(guess[3:length(guess)]))) {
      stop(paste("length sd is different of 1 or ",length(par[3:length(guess)]),sep=""))
    }
    if (prior.sd <= 0)
      stop("prior.sd must be a positive number")
  }

  ## call the Metropolis algorithm for MCMC
  chain <- try(evalWithTimeout(metrop(obj=.logpost,initial=guess,time=time,delta=delta,dist=dist,x=x,
                                      mean=prior.mean,sd=prior.sd,
                                      nbatch=(size-1)*jump+burn,blen=1,nspac=1,scale=scale),
                               timeout=max.time*60,onTimeout="error"),silent=TRUE)
  if(length(class(chain)) == 1 && class(chain) == "try-error") {
    stop("Time limit exceeded")
  }
  out$post <- chain$batch[seq(burn,length(chain$batch[,1]),jump),]

  # evaluating the point estimates
  if (x.k != 1) {
    out$par    <- c(mean(out$post[,1]),mean(out$post[,2]),apply(out$post[,3:length(out$post[1,])],2,mean))
    out$par.sd <- c(sd(out$post[,1]),    sd(out$post[,2]),apply(out$post[,3:length(out$post[1,])],2,sd))
  } else { 
    out$par    <- c(mean(out$post[,1]),mean(out$post[,2]),mean(out$post[,3:length(out$post[1,])]))
    out$par.sd <- c(sd(out$post[,1]),    sd(out$post[,2]),  sd(out$post[,3:length(out$post[1,])]))
  }
  # evaluating the interval estimates
  out$par.HPD   <- cbind(c(HPDinterval(as.mcmc(out$post[,1]),0.95)),
                         c(HPDinterval(as.mcmc(out$post[,2]),0.95)))
  for (i in 3:length(out$post[1,])) {
    out$par.HPD <- cbind(out$par.HPD,
                         c(HPDinterval(as.mcmc(out$post[,i]),0.95)))
  }

  # evaluating DIC
  aux.loglik  <- rep(0,size)
  for (j in 1:size)
    aux.loglik[j] <- c(.lik.tbs(out$post[j,],time,delta,dist,x))
  loglik <- mean(-2*aux.loglik)
  rm(aux.loglik)
  out$DIC <- 2*loglik+2*.lik.tbs(out$par,time,delta,dist,x)

  # evaluating error of the model
  aux.error    <- matrix(0,length(time),size)
  for (j in 1:size) {
      if (x.k != 1) {
        aux.error[,j]    <- (.g.lambda(log(time),out$post[j,1])-
           .g.lambda(x%*%matrix(out$post[j,3:length(out$post[1,])],length(out$post[j,3:length(out$post[1,])]),1),out$post[j,1]))
      } else {
        aux.error[,j]    <- (.g.lambda(log(time),out$post[j,1])-
                             .g.lambda(out$post[j,3]*x,out$post[j,1]))
      }
  }
  error <- rep(0,length(time))
  for (i in 1:length(time))
    error[i] <- mean(aux.error[i,])
  rm(aux.error)
  out$error <- error

  # time spent to do BE
  out$run.time <- .gettime() - initial.time
  return(out)
}

