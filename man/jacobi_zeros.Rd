\name{jacobiZeros}
\alias{jacobiZeros}
\title{Zeros of Jacobi Polynomials}

\description{
 Calculates the zeros of nth order Jacobi polynomials}
} 

\usage{
z = jacobiZeros(Q, a=0.0, b=0.0)
}
\arguments{
  \item{Q}{Number of zeros, order = Q - 1}
  \item{a}{\eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ A vector containing the zeros of the Jacobi Polynomial}


\details{
  An iterative procedure using the Zeros of Chebychev polynomials (whose
  zeros are known analitically) is used to compute the zeros
}

 



\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{

## Calculate the zeros of the 5th order Jacobi polynomials
z = jacobiZeros(5, 1, 0.5)


}

\keyword{math}



  
