\name{rcppcomdel}
\alias{rcppcomdel}
\title{Matrix computations for the full conditional distribution of the skewness nodes}
\description{
  This function performs matrix computations which are required for simulation from the full conditional distribution of the skewness latent variables, \eqn{\delta}'s.
}
\usage{
rcppcomdel(W, U, H)
}
\arguments{
  \item{W}{matrix of the latent variables \eqn{W}.}
  \item{U}{matrix of the latent variables \eqn{U}.}
  \item{H}{inverse of the spatial correlation matrix.}
}
\details{
Given the values of \eqn{W}, \eqn{U} and \eqn{H}, this function carries out the following R code in C++
\deqn{for(j in 1:n){A <- diag(U[,j]) %*% H; B <- B + c(A %*% W[,j]); C <- C + A %*% diag(U[,j])
}
}
and returns the matrices B and C.
}
\value{
  'rcppcomdel' returns the matrices which are used in the full conditional distribution of \eqn{\delta}'s.
}
\author{
  Behzad Mahmoudian \email{b.mahmoudian@qom.ac.ir}
}
