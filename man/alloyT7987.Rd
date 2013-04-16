% File TBSSurvival/man/alloyT7987.Rd
% Part of the TBSSurvival package for R (http://www.R-project.org)
% Copyright (C) 2012, Adriano Polpo, Cassio de Campos, Debajyoti Sinha
%                     Jianchang Lin and Stuart Lipsitz.
% Distributed under GPL 3 or later

\name{alloyT7987}

%%%
\alias{alloyT7987}
\title{Data set Alloy T7987}
\description{
Alloy - T7987: data extracted from Meeker and Escobar (1998), pp. 131.
}
\value{
  The two data variables are:
  \item{time}{Time of failure of the specimen.}
  \item{delta}{censoring indication (0 means right-censored, 1 means no censoring).}
}
\references{
  Meeker, W. and Escobar, L. (1998) \emph{Statistical Methods for Reliability Data}. Willey, ISBN 0-471-14328-6.
}
\seealso{
\code{\link{tbs.survreg.mle}}
}
\examples{
## See \link{tbs.survreg.mle} and \link{tbs.survreg.be}.

## The data looks like this:
## time delta
##  94     1
##  96     1
##  99     1
##  99     1
## 104     1
## 108     1
## 112     1
## 114     1
## 117     1
## 117     1
}
\keyword{TBS}
