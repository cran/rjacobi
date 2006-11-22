\name{lagrange}
\alias{lagrangeGJ}
\alias{lagrangeGLJ}
\alias{lagrangeGRJM}
\alias{lagrangeGRJP}

\title{Lagrange Interpolants of Gauss-Jacobi Quadrature}

\description{
 Calculates the ith Lagrange interpolation polinomial for  Q point Gauss-Jacobi quadrature.
} 

\usage{
y = lagrangeGJ(x, i, z, a=0, b=0)
y = lagrangeGLJ(x, i, z, a=0, b=0)
y = lagrangeGRJM(x, i, z, a=0, b=0)
y = lagrangeGRJP(x, i, z, a=0, b=0)
}
\arguments{
  \item{x}{Points where to calculate the Lagrange polynomial}
  \item{i}{ Which interpolant to calculate}
  \item{z}{ Quadrature nodes calculated by \code{\link{zerosGJ}} family
  of functions}
  \item{a}{ \eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{ \eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ The value of the Lagrangean interpolator}

\details{

  This function calculates

  \deqn{ y = h_i(x) }
  
  where x is any point in the domain -1 < x < 1 and i is the
  corresponding Lagrangean interpolator with 1 <= i <= Q.
  

      

  Different types of quadrature are possible:

  \itemize{
    \item Gauss-Jacobi, suffix \code{GJ}
    \item Gauss-Lobatto-Jacobi, suffix \code{GLJ}, includes both end
      points
    \item Gauss-Radau-Jacobi, suffic \code{GRJM}, includes the -1 end
      point
    \item Gauss-Radau-Jacobi, suffic \code{GRJP}, includes the +1 end
      point}


    This function is used internally (in C code) to calculate the
    interpolation matrix.
    
    

  
}
  

\seealso{
  \code{\link{quadrature}}
  \code{\link{zerosGJ}}
  \code{\link{interpmatGJ}}}


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{

## Plots the Lagrange interpolator for GLJ Q=6

z <- zerosGLJ(6)

x <- seq(-1, 1, len=101)

plot(x, lagrangeGLJ(x, 1, z), ty='l', xlab="x", ylab="y",
main="Lagrange polynomials, GLJ, Q=6")
lines(x, lagrangeGLJ(x, 2, z))
lines(x, lagrangeGLJ(x, 3, z))
lines(x, lagrangeGLJ(x, 4, z))
lines(x, lagrangeGLJ(x, 5, z))
lines(x, lagrangeGLJ(x, 6, z))

abline(v=z, col='red')
abline(h=1, col='red')


}

\keyword{math}
\keyword{dplot}


