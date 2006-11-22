\name{zerosGJ}
\alias{zerosGJ}
\alias{zerosGLJ}
\alias{zerosGRJM}
\alias{zerosGRJP}

\title{Nodes of Gauss-Jacobi Quadrature}

\description{
 Calculates the nodes of the Q point Gauss-Jacobi quadrature}
} 

\usage{
z = zerosGJ(Q, a=0.0, b=0.0)
z = zerosGLJ(Q, a=0.0, b=0.0)
z = zerosGRJM(Q, a=0.0, b=0.0)
z = zerosGRJP(Q, a=0.0, b=0.0)
}
\arguments{
  \item{Q}{ Number of nodes}
  \item{a}{ \eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{ \eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ A vector containing the nodes of the quadrature}


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
    

  
}

 

\seealso{
  \code{\link{quadrature}}
  \code{\link{weightsGJ}}
  \code{\link{diffmatGJ}}
  \code{\link{jacobiZeros}}}


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{

## Computes the quadrature nodes
z = zerosGJ(5, 1, 0.5)


}

\keyword{math}




  
