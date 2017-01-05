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

## This code is used for testing purposes. The TBSSurvival library does not
## depend on it for any of its functionalities

## installpacks <- function(loc=NULL,repos="http://stat.ethz.ch/CRAN/") {
##   ## set the repository to use
##   options(repos=repos)
##   ## install the packages
##   install.packages("coda",lib=loc)
##   install.packages("mcmc",lib=loc)
##   install.packages("normalp",lib=loc)
##   install.packages("R.methodsS3",lib=loc)
##   install.packages("R.oo",lib=loc)
##   install.packages("R.utils",lib=loc)
##   install.packages("Rsolnp",lib=loc)
##   install.packages("survival",lib=loc)
## #  install.packages("e1071",lib=loc)
## #  install.packages("eha",lib=loc)
##   install.packages("truncnorm",lib=loc)
##   install.packages("BMS",lib=loc)
  
##   ## this following line install the TBS package itself, so nothing else is needed.
##   ## For testing, sometimes it is better to work without installing it for a while...
##   ##      install.packages('./TBSSurvival_version.tar.gz',repos=NULL,type="source")
## }

## loadlibs <- function(libdir=NULL) {
##   w <- options("warn")
##   options("warn" = -1)
##   if (require("TBSSurvival",quietly=TRUE)==FALSE) {
##     library("BMS",lib.loc=libdir)
##     library("coda",lib.loc=libdir)
##     library("mcmc",lib.loc=libdir)
##     library("normalp",lib.loc=libdir)
##     library("R.methodsS3",lib.loc=libdir)
##     library("R.oo",lib.loc=libdir)
##     library("R.utils",lib.loc=libdir)
##     library("Rsolnp",lib.loc=libdir)
##     library("survival",lib.loc=libdir)
## #    library("e1071",lib.loc=libdir)
## #    library("eha",lib.loc=libdir)
##     library("truncnorm",lib.loc=libdir)
##     source("../R/tbs.survreg.be.r")
##     source("../R/ptbs.r")
##     source("../R/qtbs.r")
##     source("../R/dtbs.r")
##     source("../R/rtbs.r")
##     source("../R/htbs.r")
##     source("../R/tbs.survreg.mle.r")
##     source("../R/local.r")
##     source("../R/dt2.r")
##     source("../R/dlogis2.r")
##     source("../R/dist.error.r")
    
##   } else {
##     library("TBSSurvival")
##   }
##   options("warn" = w[[1]])
## }

## ## Load data
## alloyT7987 <- read.table("../data/alloyT7987.txt",header=TRUE)

library("TBSSurvival")

## loadlibs()
set.seed(1)

####################
## simple test with the GBSG2 (German Breast Cancer Group 2) data set from the ipred package
## SLOW ON CRAN
## library(ipred)
## data(GBSG2)
## cat('Running MLE on GBSG2 (from ipred package) without covariates\n')
## s=tbs.survreg.mle(survival::Surv(GBSG2$time,GBSG2$cens==1) ~ 1,verbose=TRUE) ## as optim method not given, it tries all methods
## plot(s)
## lines(s$cauchy,col=2)
## lines(s$t,col=3)
## lines(s$doubexp,col=4)
## lines(s$logistic,col=5)
## lines(s$logistic,col=5,lwd=4)

## ####################
## test with the colon data set from the survival package
## SLOW ON CRAN
library(survival)
data(colon)
cat('Running MLE on colon (from survival package) without covariates\n')
s=tbs.survreg.mle(survival::Surv(colon$time,colon$status==1) ~ 1,dist="norm",method=c("BFGS","Nelder-Mead"),verbose=TRUE,gradient=FALSE)
## plot(s)

## b=tbs.survreg.be(survival::Surv(colon$time,colon$status==1) ~ 1,dist=dist.error("norm"),burn=10000,jump=200,scale=0.05)
## b2=tbs.survreg.be(survival::Surv(colon$time,colon$status==1) ~ 1,dist=dist.error("doubexp"),burn=10000,jump=200,scale=0.05)
## plot(b)
## lines(b2,col=c(2,2))

## ## with covariate
## cat('Running MLE on colon (from survival package) with covariate=age60\n')
## colon$age60=as.numeric(colon$age>60) #threshold defined from medical papers
## s=tbs.survreg.mle(survival::Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",method=c("BFGS","Nelder-Mead"),verbose=TRUE)
## b=tbs.survreg.be(survival::Surv(colon$time,colon$status==1) ~ colon$age60,dist="norm",burn=10000,jump=200,scale=0.05)

## ## simple test with the Alloy T7987 data set (available in TBSSurvival)
## ## FAILS ON CRAN
## data(alloyT7987)

## ## MLE Estimation
## method <- c("Rsolnp","BFGS","Nelder-Mead","CG","SANN")
## for (j in 1:length(method)) {
##   for (i in 1:5) {
##     tbs.mle.logistic <- tbs.survreg.mle(survival::Surv(alloyT7987$time,alloyT7987$delta) ~ 1,dist=dist.error("logistic"),method=method[j])
##     cat("i: ",i,"  - method: ",method[j],"  - log.lik: ",tbs.mle.logistic$log.lik,"\n")
##   }
## }

