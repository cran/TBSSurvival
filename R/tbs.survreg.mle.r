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

## Maximum likelihood estimation for TBS
## By default, try all optimization methods listed below and keep the best solution
## max.time is the time limit in minutes for each method (<= 0 means no limit), ncore the number of cores
##   to use in case multicore is available, nstart the number of feasible starting
##   points to use, dist has to be one of the available distributions, currently
##   one of "norm", "t", "cauchy", "doubexp", "logistic"
## formula is a R formula with a Surv object on the left side
tbs.survreg.mle <- function(formula,dist="norm",method=c("BFGS","Rsolnp","Nelder-Mead","CG","SANN"),verbose=FALSE,nstart=10,max.time=-1,ncore=1) {
  ## this meta-method only records the elapsed time and call the max likelihood estimation function
  ## for each of the methods given until one of them converges. It is supposed that at least one method
  ## is given, and that dist is one of those implemented by tbs.survreg.
  initial.time <- .gettime()
  bestout <- NULL
  for(i in 1:length(method)) {
    ## call the estimation function
    out <- .tbs.survreg(formula,dist=dist,method=method[i],verbose=verbose,ncore=ncore,max.time=max.time,nstart=nstart)
    ## if converged, we are happy
    if(out$convergence) {
      if(is.null(bestout) || out$log.lik > bestout$log.lik) {
        bestout <- out
      }
    }
  }
  ## check if at least one method found a solution...
  if (is.null(bestout)) {
    if(verbose) cat('No method has found a solution\n')
    bestout$convergence <- FALSE
    bestout$method <- NULL
  }
  ## record the call arguments and the formula
  bestout$call <- match.call()
  bestout$formula <- formula
  ## run.time is returned in minutes
  bestout$run.time <- .gettime() - initial.time
  return(bestout)
}
