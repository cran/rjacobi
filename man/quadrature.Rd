\name{quadrature}
\alias{quadrature}

\title{Quadrature information}

\description{Builds a data structure with every quadrature related stuff
  such as integrals, derivatives and interpolations}

\usage{
quad <- quadrature(Q, type=c('GJ', 'GLJ', 'GRJM', 'GRJP', a=0, b=0,
xinterp=NULL))
}

\arguments{
  \item{Q}{Number of quadrature points}
  \item{type}{A string specifying the quadrature type}
  \item{a}{\eqn{\alpha} weight}
  \item{b}{\eqn{\beta} weight}
  \item{xinterp}{Points where the interpolation matrix should be
    calculated}
}
\value{
  \item{x}{Quadrature zeros}
  \item{a}{\eqn{\alpha} weight}
  \item{b}{\eqn{\beta} weight}
  \item{w}{Quadrature weights}
  \item{D}{Derivative matrix}
  \item{type}{Quadrature type}
  \item{lagrange}{Lagrange interpolants}
  \item{xinterp}{Interpolation points}
  \item{imat}{Interpolation matrix at \code{xinterp}}
}

\details{
  This function is a wrapper around the various functions
  \itemize{
    \item \code{zerosXXX}
    \item \code{weightsXXX}
    \item \code{diffmatXXX}
    \item \code{interpmatXXX}
    \item \code{lagrangeXXX}
  }
  
  It returns a list of class \code{quadrature} with all the above
  fields.

  The interpolation points can be modified with the
  \code{\link{xinterp}} function.
}

\seealso{
  \code{\link{lagrangeGJ}}
  \code{\link{zerosGJ}}
  \code{\link{weightsGJ}}
  \code{\link{diffmatGJ}}
  \code{\link{xinterp}}
  \code{\link{interpmatGJ}}}


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}


\examples{

## Calculates the numerical derivative of jacobi(x,4)
q = quadrature(5)
f <- jacobi(q$x, 4)
df <- q$D \%*\% f

x <- seq(-1, 1, len=51)
ff <- djacobi(x, 4)

plot(x, ff, ty='l', xlab="x", ylab="dP_4^(0,0)/dx", main="Analytical and
Numerical derivative")
points(q$x,df)


}


\keyword{math}
\keyword{list}
\keyword{array}
    
