% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convolution.R
\name{ConvDecodeHard}
\alias{ConvDecodeHard}
\title{Convolutional decoding of a code (hard decision).}
\usage{
ConvDecodeHard(code, conv.encoder = NULL, terminate = TRUE,
  punctuation.matrix = NULL, visualize = FALSE)
}
\arguments{
\item{code}{The code to be decoded.}

\item{conv.encoder}{Convolutional encoder used for encoding.}

\item{terminate}{flag If the code is terminated.}

\item{punctuation.matrix}{If not null the code is depunctured prior to the decode algorithm.}

\item{visualize}{If TRUE a beamer PDF file is generated showing the decode process.}
}
\value{
The hard-decoded message vector.
}
\description{
Decodes a convolutional codeword using hard decision decoding.
}
\examples{
coder <- ConvGenerateEncoder(2,2,c(7,5))
plain <- c(1,0,0,1,1)

# standard encoding and decoding with termination
coded <- ConvEncode(plain, coder)
ConvDecodeHard(coded, coder)

# without termination
coded <- ConvEncode(plain, coder, terminate = FALSE)
ConvDecodeHard(coded, coder, terminate = FALSE)

# with punctuation and visualization
# punctuation <- ConvGetPunctuationMatrix(c(1,1,0,1,1,0), coder)
# coded <- ConvEncode(plain, coder, punctuation.matrix = punctuation)
# ConvDecodeHard(coded, coder, punctuation.matrix = punctuation, visualize = TRUE)

# with message distortion
coded <- ConvEncode(plain, coder)
noisy <- ApplyNoise(coded, SNR.db = 3)
ConvDecodeHard(noisy, coder)
}
\author{
Martin Nocker
}

