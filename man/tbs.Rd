% File TBSSurvival/man/tbs.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{tbs}

%%%
\alias{dtbs}
\alias{ptbs}
\alias{qtbs}
\alias{rtbs}
\alias{htbs}

\title{The TBS Time Failure Distribution}
\description{
  Density, distribution function, quantile function, random
  generation and hazard function for the TBS distribution.
}
\usage{
 dtbs(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm")
 ptbs(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm")
 qtbs(p,lambda=1,xi=1,beta=1,x=NULL,dist="norm")
 rtbs(n,lambda=1,xi=1,beta=1,x=NULL,dist="norm")
 htbs(time,lambda=1,xi=1,beta=1,x=NULL,dist="norm")
}
\arguments{
  \item{time}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{lambda}{parameter of TBS.}
  \item{xi}{parameter of the error distribution.}
  \item{beta}{parameter of the linear regressor.}
  \item{x}{vector/matrix of co-variables, x=NULL if there are not co-variables.}
  \item{dist}{Distribution of error, dist = "norm", "t", "doubexp", "cauchy" or "logistic".}
}
\details{
  The density, distribution function, quantile function, random
  generation and hazard function for the failure time of a TBS Model.
  The distribution of error can be chosen from Normal, t-Student, Cauchy, Logistic and Doub-Exponential (Laplace).
}
\value{
     `dtbs' gives the density, `ptbs' gives the distribution
     function, `qtbs' gives the quantile function, `rtbs'
     generates random deviates, `htbs' gives the hazard function.
}
\examples{
  ptbs(1,lambda=2,xi=1,beta=1,dist="norm")
}
