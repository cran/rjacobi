\name{quadInterp}
\alias{quadInterp}
\title{Polynomial Interpolation}

\description{
  Interpolates a function given at the quadrature nodes
} 

\usage{
df = quadInterp(imat, f)
}
\arguments{
  \item{imat} {An interpolation matrix created by \code{\link{interpmatGJ}} family of functions}
  \item{f} {A function or a vector with function values at quadrature points}
}
\value{ Returns the interpolated values of f}



\seealso{
  \code{\link{quadrature}}
  \code{\link{xinterp}}
  \code{\link{zerosGJ}}
  \code{\link{interpmatGJ}}
  \code{\link{quadGJ}}}


  



\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}


\keyword{math}



  
