\name{xinterp}
\alias{xinterp}
\alias{xinterp<-}
\title{Quadrature interpolation}

\description{Returns and modifies a interpolation points of a quadrature}

\usage{
x <- xinterp(q)
xinterp(q) <- value
}

\arguments{
  \item{q}{An object of class quadrature returned by the quadrature
    function}
  \item{value}{Points where the interpolation should be calculated}
}
\value{
  Interpolation points
}

\details{
  This function is used to set or modify the points where interpolation
  occurs. It modifies the fields \code{xinterp} and \code{imat} of the
  an object with class \code{quadrature}.
}

\seealso{
  \code{\link{lagrangeGJ}}
  \code{\link{quadrature}}
  \code{\link{quadInterp}}
  \code{\link{interpmatGJ}}}


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}


\examples{

## Interpolates the Runge function in the interval (-1,1)

runge <- function(x){1 / (1 + 25*x^2)}

q <- quadrature(13, type='GLJ')


z <- q$x

x <- seq(-1, 1, by=0.02)
xinterp(q) <- x

x2 <- seq(-1, 1, len=201)

fe <- runge(x2)

f <- runge(z)

finterp <- quadInterp(q, runge)


plot(x2, fe, ty='l', xlab='x', ylab='y',
main='Runge function interpolation')

points(z, f)
lines(x, finterp, lty=2)



}
\keyword{math}
\keyword{list}

    
