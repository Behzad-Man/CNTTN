\name{rcppng}
\alias{rcppng}
\title{Simulation from the full conditional distribution of the spatial random effect}
\description{
  This function performs matrix computations for the full conditional distribution of \eqn{W}.
}
\usage{
rcppng(sigma, sigw, isigt, invL, muW, rngN)
 }
\arguments{
  \item{sigma}{inverse of the spatial correlation matrix.}
  \item{sigw}{value of the scale parameter \eqn{\gamma}.}
  \item{isigt}{inverse of the error variance \eqn{\tau^2}.}
  \item{invL}{vector of inverse of the \eqn{\lambda}'s values.}
  \item{muW}{mean part of the full conditional distribution.}
  \item{rngN}{matrix of N(0,1) generated random numbers.}
}
\details{
Given the arguments values, this function implements the following R code in C++

	\deqn{for(j in 1:n){ K.W <- chol2inv(chol(sigma + sigw*sigw*diag(invL[,j])*isigt));
	W[,j] <- K.W %*% muW[,j] + c(rngN %*% chol(K.W))
	}
}
}
\value{
  'rcppng' returns the simulated latent variables \eqn{W}'s from theirs full conditional distributions.
}
\author{
  Behzad Mahmoudian \email{b.mahmoudian@qom.ac.ir}
}
