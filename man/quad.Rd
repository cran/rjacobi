\name{quad}
\alias{quadGJ}
\alias{quadGLJ}
\alias{quadGRJM}
\alias{quadGRJP}

\title{Gauss-Jacobi Quadrature Parameters}

\description{
 Calculates the Interpolation matrix of the Q point Gauss-Jacobi quadrature.
} 

\usage{
quad <- quadGJ(Q, a=0, b=0)
quad <- quadGLJ(Q, a=0, b=0)
quad <- quadGRJM(Q, a=0, b=0)
quad <- quadGRJP(Q, a=0, b=0)
}
\arguments{
  \item{Q}{Number of quadrature points}
  \item{a}{\eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{
  \item{x}{Quadrature nodes}
  \item{w}{Quadrature weights}
  \item{D}{Derivative matrix}
  \item{interp}{A function that computes the interpolation matrix given
  a set of nodes.}
}
\details{

  This function is actually a set of wrappers around the functions
  \code{\link{zerosGJ}}, \code{\link{weightsGJ}},
  \code{\link{diffmatGJ}} and \code{\link{interpmatGJ}}.

  It is basically a friendly interface to these functions.
      

  Different types of quadrature are possible:

  \itemize{
    \item Gauss-Jacobi, suffix \code{GJ}
    \item Gauss-Lobatto-Jacobi, suffix \code{GLJ}, includes both end
      points
    \item Gauss-Radau-Jacobi, suffic \code{GRJM}, includes the -1 end
      point
    \item Gauss-Radau-Jacobi, suffic \code{GRJP}, includes the +1 end
      point}


    
    

  
}
  

\seealso{
  \code{\link{quadrature}}
  \code{\link{zerosGJ}}
  \code{\link{weightsGJ}}
  \code{\link{diffmatGJ}}
  \code{\link{interpmatGJ}}}
  


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}


\examples{

## Calculates the numerical derivative of jacobi(x,4)
q = quadGJ(5)
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


