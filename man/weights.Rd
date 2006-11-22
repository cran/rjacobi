\name{weightsGJ}
\alias{weightsGJ}
\alias{weightsGLJ}
\alias{weightsGRJM}
\alias{weightsGRJP}

\title{Weights of Gauss-Jacobi Quadrature}

\description{
 Calculates the weights of the Q point Gauss-Jacobi quadrature.
} 

\usage{
z = weightsGJ(z, a=0.0, b=0.0)
z = weightsGLJ(z, a=0.0, b=0.0)
z = weightsGRJM(z, a=0.0, b=0.0)
z = weightsGRJP(z, a=0.0, b=0.0)

}
\arguments{
  \item{z}{Quadrature node calculated with \code{\link{zerosGJ}}
  according to quadrature type}
  \item{a}{\eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ A vector containing the weights of the quadrature}


\details{

  Different types of quadrature are possible:

  \itemize{
    \item Gauss-Jacobi, suffix \code{GJ}
    \item Gauss-Lobatto-Jacobi, suffix \code{GLJ}, includes both end
      points
    \item Gauss-Radau-Jacobi, suffic \code{GRJM}, includes the -1 end
      point
    \item Gauss-Radau-Jacobi, suffic \code{GRJP}, includes the +1 end
      point}

  These functions find the quadrature nodes that will be used in
  subsequent calculations.

  With these weights, integrals can be calculated according to the
  following expression:

  \deqn{\int_{-1}^{1} f(x)dx \approx \sum_{i=1}{Q}w_i f(x_i)}
  

  
}
  
\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{
### This example will shoud the exponential convergence when integrating
### smooth functions.

Q <- 2:15
f <- function(x)cos(5*x)
n <- length(Q)
integr <- double(n)

for (i in 1:n){
  z <- zerosGJ(Q[i])
  w <- weightsGJ(z)
  integr[i] <- sum(w * f(z))
}
err <- abs(integr - 2*sin(5)/5)

plot(Q, err, ty='b', log='y', main="Error of integral of cos(5*x)",
xlab="No Quadrature points", ylab="Error")




}
\seealso{
  \code{\link{zerosGJ}}
  \code{\link{diffmatGJ}}}
  
  
\keyword{math}



  
