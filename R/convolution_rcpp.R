# convolution_rcpp.R
#
# interface for convolutional codes
# provided functions:
#  - generate nsc coder
#  - generate rsc coder
#  - encode
#  - decode (soft in & out)
#  - decode (hard decision)

#library(Rcpp)
#sourceCpp('src/convolution.cpp')

#' generate nsc encoder
#'
#' generates a convolutional encoder for nonsystematic convolutional codes (nsc)
#' @author Martin Nocker
#' @param N numer ob output symbols per input symbol
#' @param M constraint length (memory length of the encoder)
#' @param generators vector of generator polynoms (one for each output symbol)
#' @return a convolutional encoder represented as a list containing:
#' N, M, 4 matrices: nextState, previousState, output, termination, nsc (flag)
#' @examples generateConvEncoder_nsc(2,2,c(7,5))
#' @export
#' @useDynLib channelcoding
#' @importFrom Rcpp sourceCpp
GenerateNscEncoder <- function(N, M, generators) {

  # nsc requires N generator polynoms
  if (length(generators) < N) {
    # stop execution if too few generators
    stop("Too few generator polynoms!")
  } else if (length(generators) > N) {
    # just use first N polynoms if too many are provided
    generators <- head(generators, N)
  }

  matrix.list <- c_generateMatrices_nsc(N,M,generators)

  conv.encoder <- list(N = N,
                       M = M,
                       generators = generators,
                       next.state = matrix.list$next.state,
                       prev.state = matrix.list$prev.state,
                       output = matrix.list$output,
                       nsc = TRUE,
                       termination = vector())

  return(conv.encoder)
}

#' generate rsc encoder
#'
#' generates a convolutional encoder for recursive systematic codes (rsc)
#' @author Martin Nocker
#' @param N numer ob output symbols per input symbol
#' @param M constraint length (memory length of the encoder)
#' @param generators vector of generator polynoms (one for each output symbol and one for the recursion)
#' @return a convolutional encoder represented as a list containing:
#' N, M, 4 matrices: nextState, previousState, output, termination, nsc (flag)
#' @export
#' @useDynLib channelcoding
#' @importFrom Rcpp sourceCpp
GenerateRscEncoder <- function(N, M, generators) {

  # rsc requires N+1 generator polynoms
  if (length(generators) < N+1) {
    # stop execution if too few generators
    stop("Too few generator polynoms!")
  } else if (length(generators) > N+1) {
    # just use first N polynoms if too many are provided
    generators <- head(generators, N+1)
  }

  matrix.list <- c_generateMatrices_rsc(N,M,generators)

  conv.encoder <- list(N = N,
                       M = M,
                       generators = generators,
                       next.state = matrix.list$next.state,
                       prev.state = matrix.list$prev.state,
                       output = matrix.list$output,
                       nsc = FALSE,
                       termination = matrix.list$termination)

  return(conv.encoder)
}

#' convolutional encoding of a message
#'
#' \code{ConvEncode} produces a convolutional code of a message based on the encoder passed as an argument
#' @author Martin Nocker
#' @param message the message to be encoded
#' @param conv.encoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code should be terminated, default: TRUE
#' @return the encoded message
#' @export
#' @useDynLib channelcoding
#' @importFrom Rcpp sourceCpp
ConvEncode <- function(message, conv.encoder, terminate = TRUE) {

  code <- c_convolutionEncode(message,
                              conv.encoder$N,
                              conv.encoder$M,
                              conv.encoder$next.state,
                              conv.encoder$output,
                              as.integer(conv.encoder$nsc),
                              conv.encoder$termination,
                              as.integer(terminate))
  return(code)
}

#' convolutional decoding of a code (soft decision)
#'
#' \code{ConvDecode} decodes a codeword that was encoded with the given encoder.
#' This decoder is a soft-input soft-output decoder
#' @author Martin Nocker
#' @param code the code to be decoded
#' @param conv.encoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code was terminated, default: TRUE
#' @return the decoded message, list(softOutput, hardOutput)
#'
#' @useDynLib channelcoding
#' @export
#' @export
#' @useDynLib channelcoding
#' @importFrom Rcpp sourceCpp
ConvDecode <- function(code, conv.encoder, terminate = TRUE) {

  output <- c_convolutionDecode(code,
                                conv.encoder$N,
                                conv.encoder$M,
                                conv.encoder$prev.state,
                                conv.encoder$output)

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    soft <- head(output$soft.output, length(output$soft.output) - conv.encoder$M)
    hard <- head(output$hard.output, length(output$hard.output) - conv.encoder$M)

    newlist <- list(soft.output = soft, hard.output = hard)
    return(newlist)
  }

  return(output)
}

#' convolutional decoding of a code (hard decision)
#'
#' \code{ConvDecodeHard} decodes a codeword that was encoded with the given encoder.
#' This decoder is a hard-decision decoder
#' @author Martin Nocker
#' @param code the code to be decoded
#' @param conv.encoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code was terminated, default: TRUE
#' @return the hard-decoded message vector
#' @export
#' @useDynLib channelcoding
#' @importFrom Rcpp sourceCpp
ConvDecodeHard <- function(code, conv.encoder, terminate = TRUE) {

  output <- c_convolutionDecode_hard(code,
                                     conv.encoder$N,
                                     conv.encoder$M,
                                     conv.encoder$prev.state,
                                     conv.encoder$output)

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    return(head(output, length(output) - conv.encoder$M))
  }

  return(output)
}
