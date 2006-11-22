\name{jacobi}
\alias{jacobi}
\title{Jacobi Polynomials}

\description{
  Computes Jacobi polynomials
} 

\usage{
y = jacobi(x,n,a=0,b=0)
}
\arguments{
  \item{x}{A vector containing the points where the polynomials  should be calculated. \eqn{-1\le x_i \le 1}}
  \item{n}{Order of Jacobi polynomial to be calculated}
  \item{a}{\eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0.}
}
\value{ A vector containing the values of the computed Jacobi
  Polynomials}


\details{
  The Jacobi Polynomials are computed from the recurrence formula:

  \deqn{a^1_n P_{n+1}^{\alpha,\beta}(x) = (a^2_n + a^3_n
    x)P_n^{\alpha,\beta}(x) - a^4_n P_{n-1}^{\alpha,\beta}(x)}
}

 



\references{
  Abramowitz, Milton and Stegun, Irene (editors); "Handbook of Mathematical functions",
  Dover Publications 1965.
  
  Karniadakis, George Em and Sherwin, Spencer; "Spectral/hp Element Methods for Computational
  Fluid Dynamics", Oxford Science Publications, 2nd edition, 2005.
}

\examples{

## Example to plot the 5th order Jacobi Polynomial
x <- seq(-1, 1, len=101)
y <- jacobi(x, 5, 0, 0)
plot(x, y, ty='l')



}

\keyword{math}
\keyword{arith}


  
