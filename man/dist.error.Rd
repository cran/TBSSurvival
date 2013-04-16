% File TBSSurvival/man/dist.error.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2013, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{dist.error}
\alias{dist.error}
%%%

\title{The Wrapper for error distribution functions to use with TBS}
\description{
  Constructor of a list of density function, distribution function, quantile function, random
  generation and hazard function for the Transform-Both-Sides (TBS) distribution.
}
\usage{
  dist.error(dist="norm")
}
\arguments{
  \item{dist}{Distribution of error, dist = "norm", "t", "doubexp", "cauchy" or "logistic". A choice "all" can also be
given, in which case a list with all the available distributions will be generated.}
}
\details{
  This wrapper creates a list with the necessary information about a distribution to be used with the Transform-Both-Sides (TBS) model.
  It returns a list with five items, namely the density function, distribution function, quantile function, random generation function,
   and name (a string). For example, with dist = "norm" it gives list(dnorm,pnorm,qnorm,rnorm,"norm"). The idea is that
  the user can implement its own list to be used with the TBS, as long as it is a zero-centered unimodal symmetric distribution, and replace the call
  of dist.error with their own list (given in a similar way as the output of dist.error).
}
\value{
     A list(d.dist,p.dist,q.dist,r.dist,name.dist) according with the chosen distribution.
     The functions d.dist,p.dist,q.dist,r.dist have to accept exactly two arguments: the first is the actual argument to the function, and the second is a parameter. 
}
\examples{
  ## this will return list(dcauchy,pcauchy,qcauchy,rcauchy,"cauchy")
  dist = dist.error("cauchy")

  ## a user-built distribution would look like:
  ## dist = list(
  ##         d = function(x,xi) dmydistrib(x,param=xi), # density
  ##         p = function(x,xi) pmydistrib(x,param=xi), # distr
  ##         q = function(x,xi) qmydistrib(x,param=xi), # quantile
  ##         r = function(x,xi) rmydistrib(x,param=xi), # generation
  ##         name = "mydistrib"
  ##        )
}
