
\name{djacobi}
\alias{djacobi}
\title{Derivative of Jacobi Polynomials}

\description{
 Calculates the derivative of the nth order Jacobi polynomials at an array of np points given by the vector x.
} 

\usage{
y = djacobi(x, n, a=0.0, b=0.0)
}
\arguments{
  \item{x}{ A vector of length np containing the points where the polynomials  should be calculated. \eqn{-1\le x_i \le 1}}
  \item{n}{Order of Jacobi polynomial to be calculated}
  \item{a}{\eqn{\alpha} parameter of Jacobi polynomial. Defaults to 0}
  \item{b}{\eqn{\beta} parameter of Jacobi polynomial. Defaults to 0}
}
\value{ A vector containing the values of the computed derivatives of Jacobi
  Polynomials}


\details{
  The derivative of the Jacobi Polynomials is computed from the
  following equation:
  \deqn{\frac{d}{dx}P^{\alpha,\beta}_n(x) = \frac{1}{2}(\alpha + \beta +
    n + 1) P_n^{\alpha+1,\beta+1}(x)}

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
y <- djacobi(x, 5, 0, 0)
plot(x, y, ty='l')

}

\keyword{math}
\keyword{arith}


