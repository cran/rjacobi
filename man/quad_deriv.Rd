\name{quadDeriv}
\alias{quadDeriv}
\title{Quadrature Derivation}

\description{
  Derivates a function given at the quadrature nodes
} 

\usage{
df = quadDeriv(quad, f)
}
\arguments{
  \item{quad}{ A list created by the \code{\link{quadGJ}} family of functions}
  \item{f}{A function or a vector with function values at quadrature points}
}
\value{ Returns the approximate derivative of the function at quadrature
nodes}


\seealso{
  \code{\link{quadrature}}
  \code{\link{zerosGJ}}
  \code{\link{diffmatGJ}}
  \code{\link{quadGJ}}}

 



\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}


\keyword{math}



  
