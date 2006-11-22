\name{interpmat}
\alias{interpmatGJ}
\alias{interpmatGLJ}
\alias{interpmatGRJM}
\alias{interpmatGRJP}

\title{Interpolation Matrix of Gauss-Jacobi Quadrature}

\description{
 Calculates the Interpolation matrix of the Q point Gauss-Jacobi quadrature.
} 

\usage{
Imat = interpmatGJ(zp, x, a=0.0, b=0.0)
Imat = interpmatGLJ(zp, x, a=0.0, b=0.0)
Imat = interpmatGRJM(zp, x, a=0.0, b=0.0)
Imat = interpmatGRJP(zp, x, a=0.0, b=0.0)
}
\arguments{
  \item{zp}{Points where to calculate the interpolation}
  \item{x}{Quadrature node calculated with \code{\link{zerosGJ}}
  (depending on quadrature type)}
  \item{a}{\eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ The interpolation matrix}

\details{

  The interpolation matrix is used to interpolate functions known at
  quadrature nodes to other points. It is actually a polynomial
  interpolation using Lagrange interpolants through quadrature
  nodes. The interpolation matrix is given by:

  \deqn{I_{ij} =  h_j(x_i)}

  where \eqn{x_i} are the interpolation points and \eqn{h_j(x)} are the
  Lagrangean interpolators associated to the quadrature points.
   

      

  Different types of quadrature are possible:

  \itemize{
    \item Gauss-Jacobi, suffix \code{GJ}
    \item Gauss-Lobatto-Jacobi, suffix \code{GLJ}, includes both end
      points
    \item Gauss-Radau-Jacobi, suffic \code{GRJM}, includes the -1 end
      point
    \item Gauss-Radau-Jacobi, suffic \code{GRJP}, includes the +1 end
      point}


  With this matrix, a function known at quadrature points can be
  calculated with

  \deqn{[f]_{interp} = \{I_{ij} \}[f]_{quad} }
  

  \code{ df = D \%*\% f(x)}
  where x is a vector containning the quadrature nodes.

    
    

  
}
  

\seealso{
  \code{\link{quadrature}}
  \code{\link{xinterp}}
  \code{\link{zerosGJ}}
  \code{\link{diffmatGJ}}
  \code{\link{weightsGJ}}}


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{

## Interpolates the Runge function in the interval (-1,1)

runge <- function(x){1 / (1 + 25*x^2)}


z <- zerosGLJ(13)
x <- seq(-1, 1, by=0.02)

x2 <- seq(-1, 1, len=201)

fe <- runge(x2)

f <- runge(z)
Imat <- interpmatGLJ(x, z)

finterp <- Imat \%*\% f


plot(x2, fe, ty='l', xlab='x', ylab='y',
main='Runge function interpolation')

points(z, f)
lines(x, finterp, lty=2)



}

\keyword{math}
\keyword{array}


