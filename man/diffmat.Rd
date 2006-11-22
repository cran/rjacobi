\name{diffmatGJ}
\alias{diffmatGJ}
\alias{diffmatGLJ}
\alias{diffmatGRJM}
\alias{diffmatGRJP}

\title{Derivative Matrix of Gauss-Jacobi Quadrature}

\description{
 Calculates the derivative matrix of the Q point Gauss-Jacobi quadrature.
} 

\usage{
D = diffmatGJ(z, a=0.0, b=0.0)
D = diffmatGLJ(z, a=0.0, b=0.0)
D = diffmatGRJM(z, a=0.0, b=0.0)
D = diffmatGRJP(z, a=0.0, b=0.0)
}
\arguments{
  \item{z}{ Quadrature node calculated with \code{\link{zerosGJ}}
  (depending on quadrature type)}
  \item{a}{ \eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0 }
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ The derivative matrix}

\details{

  The derivative matrix is used to calculate derivatives at nodal points
  of functions known at the same nodal points. It is defined according
  to the following expression:
  

  \deqn{ \left(\frac{du(x)}{dx}\right)_{x=x_i} = \sum_{j=0}^{Q-1} D_{ij}
    u(x_j)}
  
  In this equation,
  \deqn{D_{ij} = \left(\frac{dh_j(x)}{dx}\right)_{x=x_i}}
  

  where \eqn{h_j(x)} is the Lagrange polynomial through the jth quadrature node
      

  Different types of quadrature are possible:

  \itemize{
    \item Gauss-Jacobi, suffix \code{GJ}
    \item Gauss-Lobatto-Jacobi, suffix \code{GLJ}, includes both end
      points
    \item Gauss-Radau-Jacobi, suffic \code{GRJM}, includes the -1 end
      point
    \item Gauss-Radau-Jacobi, suffic \code{GRJP}, includes the +1 end
      point}


  With the derivative matrix, the derivative at quadrature points can be calculated:

  \code{ df = D \%*\% f(x)}
  where x is a vector containning the quadrature nodes.

    
    

  
}
  

\seealso{
  \code{\link{zerosGJ}}
  \code{\link{weightsGJ}}
  \code{\link{quadrature}}}


\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{

## Calculates the numerical derivative of jacobi(x,4)
z = zerosGJ(5)
D = diffmatGJ(z)
f <- jacobi(z, 4)
df <- D \%*\% f

x <- seq(-1, 1, len=51)
ff <- djacobi(x, 4)

plot(x, ff, ty='l', xlab="x", ylab="dP_4^(0,0)/dx", main="Analytical and
Numerical derivative")
points(z,df)


}

\keyword{math}
\keyword{array}

