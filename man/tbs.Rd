% File TBSSurvival/man/tbs.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012-2013, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
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
  Density function, distribution function, quantile function, random
  generation function and hazard function for the Transform-Both-Sides (TBS) model.
}
\usage{
 dtbs(time,lambda=1,xi=1,beta=1,x=NULL,dist=dist.error("norm"))
 ptbs(time,lambda=1,xi=1,beta=1,x=NULL,dist=dist.error("norm"))
 qtbs(p,lambda=1,xi=1,beta=1,x=NULL,dist=dist.error("norm"))
 rtbs(n,lambda=1,xi=1,beta=1,x=NULL,dist=dist.error("norm"))
 htbs(time,lambda=1,xi=1,beta=1,x=NULL,dist=dist.error("norm"))
}
\arguments{
  \item{time}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{lambda}{parameter of TBS.}
  \item{xi}{parameter of the error distribution.}
  \item{beta}{parameter of the linear regressor.}
  \item{x}{vector/matrix of co-variables, x=NULL if there are not co-variables.}
  \item{dist}{Distribution of error, it can be string such as dist = "norm", "t", "doubexp", "cauchy" or "logistic, or it can also be
              given as a list of functions (density, distribution, quantile, random generation, name). Details below.}
}
\details{
  The density function, distribution function, quantile function, random
  generation and hazard function for the failure time of a TBS Model.
  The distribution of error can be chosen from Normal, t-Student, Cauchy, Logistic and Doub-Exponential (Laplace), or can
  be given by the user (as long as it is zero-centered, unimodal and symmetric -- TBS does not check it). See the help of \link{dist.error}
  for examples.
}
\value{
     `dtbs' gives the density, `ptbs' gives the distribution
     function, `qtbs' gives the quantile function, `rtbs'
     generates random deviates, `htbs' gives the hazard function.
}
\seealso{
\code{\link{dist.error}}
}
\examples{
  ptbs(1,lambda=2,xi=1,beta=1,dist=dist.error("norm"))
}
