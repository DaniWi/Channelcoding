% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blockcodes.R
\name{BlockSimulation}
\alias{BlockSimulation}
\title{Block Simulation.}
\usage{
BlockSimulation(coder = NULL, msg.length = 100, min.db = 0.1,
  max.db = 2, db.interval = 0.1, iterations.per.db = 100,
  visualize = FALSE)
}
\arguments{
\item{coder}{Block coder used for the simulation. Can be created via
\code{\link{BlockGenerateEncoderBCH}} or \code{\link{BlockGenerateEncoderHamming}}.}

\item{msg.length}{Message length of the randomly created messages to be encoded. Will be shortened to a multiple of the coders blocksize}

\item{min.db}{Minimum SNR to be tested.}

\item{max.db}{Maximum SNR to be tested.}

\item{db.interval}{Step between two SNRs tested.}

\item{iterations.per.db}{Number of encode and decode processes per SNR.}

\item{visualize}{If true a PDF report is generated.}
}
\value{
Dataframe containing the bit-error-rates for each SNR tested.
}
\description{
Simulation of a block encode and decode process over a noisy channel.
}
\examples{
# use all default parameters
BlockSimulation()

# Custom coder
coder <- BlockGenerateEncoderHamming(15,11)
BlockSimulation(coder, 15, 0.01, 1, 0.05, 50, FALSE)
}
\author{
Benedikt Wimmer
}

