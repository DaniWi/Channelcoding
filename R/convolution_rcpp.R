# convolution_rcpp.R
#
# interface for convolutional codes
# provided functions:
#  - generate nsc coder
#  - generate rsc coder
#  - encode
#  - decode (soft in & out)
#  - decode (hard decision)

library(Rcpp);
sourceCpp('src/convolution.cpp');

#' generate nsc encoder
#'
#' generates a convolutional encoder for nonsystematic convolutional codes (nsc)
#' @author Martin Nocker
#' @param N numer ob output symbols per input symbol
#' @param M constraint length (memory length of the encoder)
#' @param generators vector of generator polynoms (one for each output symbol)
#' @return a convolutional encoder represented as a list containing:
#' N, M, 4 matrices: nextState, previousState, output, termination, nsc (flag)
#' @example generateConvEncoder_nsc(2,2,c(7,5))
generateConvEncoder_nsc <- function(N, M, generators) {

   # nsc requires N generator polynoms
   if (length(generators) < N) {
      # stop execution if too few generators
      stop("Too few generator polynoms!");
   }
   else if (length(generators) > N) {
      # just use first N polynoms if too many are provided
      generators <- head(generators, N);
   }

   matrixList <- c_generateMatrices_nsc(N,M,generators);

   convEncoder <- list(N = N,
                       M = M,
                       generators = generators,
                       nextState = matrixList$nextState,
                       prevState = matrixList$prevState,
                       output = matrixList$output,
                       nsc = TRUE,
                       termination = vector());

   return(convEncoder);
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
#' @example generateConvEncoder_rsc(2,2,c(1,10,13))
generateConvEncoder_rsc <- function(N, M, generators) {

   # rsc requires N+1 generator polynoms
   if (length(generators) < N+1) {
      # stop execution if too few generators
      stop("Too few generator polynoms!");
   }
   else if (length(generators) > N+1) {
      # just use first N polynoms if too many are provided
      generators <- head(generators, N+1);
   }

   matrixList <- c_generateMatrices_rsc(N,M,generators);

   convEncoder_rsc <- list(N = N,
                           M = M,
                           generators = generators,
                           nextState = matrixList$nextState,
                           prevState = matrixList$prevState,
                           output = matrixList$output,
                           nsc = FALSE,
                           termination = matrixList$termination);

   return(convEncoder_rsc);
}

#' convolutional encoding of a message
#'
#' \code{conv_encode} produces a convolutional code of a message based on the encoder passed as an argument
#' @author Martin Nocker
#' @param message the message to be encoded
#' @param convEncoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code should be terminated, default: TRUE
#' @return the encoded message
#' @example
#' encoder <- generateConvEncoder_nsc(2,2,c(7,5))
#' code <- conv_encode(c(1,0,1), encoder)
conv_encode <- function(message, convEncoder, terminate = TRUE) {

   code <- c_convolutionEncode(message,
                               convEncoder$N,
                               convEncoder$M,
                               convEncoder$nextState,
                               convEncoder$output,
                               as.integer(convEncoder$nsc),
                               convEncoder$termination,
                               as.integer(terminate));
   return(code);
}

#' convolutional decoding of a code (soft decision)
#'
#' \code{conv_decode} decodes a codeword that was encoded with the given encoder.
#' This decoder is a soft-input soft-output decoder
#' @author Martin Nocker
#' @param code the code to be decoded
#' @param convEncoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code was terminated, default: TRUE
#' @return the decoded message, list(softOutput, hardOutput)
#' @example
#' encoder <- generateConvEncoder_nsc(2,2,c(7,5))
#' code <- conv_encode(c(1,0,1), encoder)
#' msg <- conv_decode(code, encoder)
conv_decode <- function(code, convEncoder, terminate = TRUE) {

   output <- c_convolutionDecode(code,
                                 convEncoder$N,
                                 convEncoder$M,
                                 convEncoder$prevState,
                                 convEncoder$output);

   # if terminated, termination bits are thrown away
   if (terminate == TRUE) {
      soft <- head(output$softOutput, length(output$softOutput) - convEncoder$M);
      hard <- head(output$hardOutput, length(output$hardOutput) - convEncoder$M);

      newlist <- list(softOutput = soft, hardOutput = hard);
      return(newlist);
   }

   return(output);
}

#' convolutional decoding of a code (hard decision)
#'
#' \code{conv_decode} decodes a codeword that was encoded with the given encoder.
#' This decoder is a hard-decision decoder
#' @author Martin Nocker
#' @param code the code to be decoded
#' @param convEncoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code was terminated, default: TRUE
#' @return the hard-decoded message vector
#' @example
#' encoder <- generateConvEncoder_nsc(2,2,c(7,5))
#' code <- conv_encode(c(1,0,1), encoder)
#' msg <- conv_decode_hard(code, encoder)
conv_decode_hard <- function(code, convEncoder, terminate = TRUE) {

   output <- c_convolutionDecode_hard(code,
                                      convEncoder$N,
                                      convEncoder$M,
                                      convEncoder$prevState,
                                      convEncoder$output);

   # if terminated, termination bits are thrown away
   if (terminate == TRUE) {
      return(head(output, length(output) - convEncoder$M));
   }

   return(output);
}
