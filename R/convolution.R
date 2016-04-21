# convolution.R
#
# interface for convolutional codes
# provided functions:
#  - generate non-recursive coder
#  - generate rsc coder
#  - encode
#  - decode (soft in & out)
#  - decode (hard decision)

#' generate convolutional encoder
#'
#' Generates a convolutional encoder for nonrecursive convolutional codes.
#' @details N is an integer and gives the number of output bits per input bit.
#'     N has to be at least two. M is an integer and gives the memory length
#'     of the encoder (number of shift register elements in the circuit). M
#'     has to be at least one. M also defines the constraint length which is
#'     M+1.
#'     The generator polynoms define how the output bits are computed for each
#'     of the N output signals. The polynoms are octal numbers. For example
#'     given a M = 2 encoder with a generator polynom of 5 for a certain output.
#'     Octal 5 means binary 101. The MSB handles the input signal, the LSB
#'     handles the output of the last memory element (last shift register
#'     element). Therefore octal 5 means the output symbol is computed as
#'     the xor combination of the input symbol and the last memory element's
#'     output.
#' @param N numer ob output symbols per input symbol
#' @param M memory length of the encoder
#' @param generators vector of N octal generator polynoms
#'     (one for each output symbol)
#' @return a convolutional encoder represented as a list containing:
#'     N, M, vector of generator polynoms,
#'     4 matrices: nextState, previousState, output and termination, rsc (flag)
#' @examples GenerateConvEncoder(2,2,c(7,5))
#' @author Martin Nocker
#' @export
GenerateConvEncoder <- function(N, M, generators) {

  stopifnot(N > 1, M > 0)

  # encoder requires N generator polynoms
  if (length(generators) < N) {
    # stop execution if too few generators
    stop("Too few generator polynoms!")
  } else if (length(generators) > N) {
    # just use first N polynoms if too many are provided
    generators <- head(generators, N)
  }

  if (!isOctal(generators)) {
    # only octal generators are accepted
    stop("At least one generator is not in octal form!")
  }

  max.generator.octal = decimalToOctal(2^(M+1) - 1)

  if (any(generators > max.generator.octal)) {
    stop("At least one generator is greater than the maximum generator!")
    # generators = maskGenerators(generators, max.generator.octal)
  }

  if (isCatastrophicEncoder(generators)) {
    warning("The result will be a catastrophic encoder!")
  }

  matrix.list <- c_generateMatrices(N,M,generators)

  conv.encoder <- list(N = N,
                       M = M,
                       generators = generators,
                       next.state = matrix.list$next.state,
                       prev.state = matrix.list$prev.state,
                       output = matrix.list$output,
                       rsc = FALSE,
                       termination = vector())

  return(conv.encoder)
}

#' generate rsc encoder
#'
#' Generates a recursive systematic convolutional (rsc) encoder.
#' @details N is an integer and gives the number of output bits per input bit.
#'     N has to be at least two. M is an integer and gives the memory length
#'     of the encoder (number of shift register elements in the circuit). M
#'     has to be at least one.
#'     The generator polynoms define how the output bits are computed for each
#'     of the N output signals. The polynoms are octal numbers. See details of
#'     \code{\link{GenerateConvEncoder}} for an example.
#'     An rsc encoder has exactly one fixed systematic output signal.
#'     The generator polynom for the systematic output doesn't have to be
#'     passed as an argument. So the generators argument contains all polynoms
#'     for non-systematic outputs and at the last position the recursion
#'     polynom. The LSB of the recursion polynom handles the input signal,
#'     the other bits handle the memory outputs. The LSB of the output polynoms
#'     handle the recursion output(!), not the original input signal. The other
#'     bits also handle the memory outputs.
#' @param N numer ob output symbols per input symbol
#' @param M memory length of the encoder
#' @param generators vector of generator polynoms
#'     (one for each non-systematic output symbol and one for the recursion)
#' @return a convolutional encoder represented as a list containing:
#'     N, M, vector of generator polynoms,
#'     4 matrices: nextState, previousState, output and termination, rsc (flag)
#' @examples GenerateRscEncoder(2,2,c(5,7))
#' @author Martin Nocker
#' @export
GenerateRscEncoder <- function(N, M, generators) {

  stopifnot(N > 1, M > 0)

  # encoder requires N generator polynoms
  if (length(generators) < N) {
    # stop execution if too few generators
    stop("Too few generator polynoms!")
  } else if (length(generators) > N) {
    # just use first N polynoms if too many are provided
    generators <- head(generators, N)
  }

  if (!isOctal(generators)) {
    # only octal generators are accepted
    stop("At least one generator is not in octal form!")
  }

  max.generator.octal = decimalToOctal(2^(M+1) - 1)

  if (any(generators > max.generator.octal)) {
    stop("At least one generator is greater than the maximum generator!")
    # generators = maskGenerators(generators, max.generator.octal)
  }

  if (isCatastrophicEncoder(generators)) {
    warning("The result will be a catastrophic encoder!")
  }

  matrix.list <- c_generateMatrices_rsc(N,M,generators)

  conv.encoder <- list(N = N,
                       M = M,
                       generators = generators,
                       next.state = matrix.list$next.state,
                       prev.state = matrix.list$prev.state,
                       output = matrix.list$output,
                       rsc = TRUE,
                       termination = matrix.list$termination)

  return(conv.encoder)
}

#' convolutional encoding of a message
#'
#' \code{ConvEncode} produces a convolutional code of a message
#' @param message the message to be encoded
#' @param conv.encoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code should be terminated, default: TRUE
#' @return the encoded message
#' @examples
#' coder <- GenerateConvEncoder(2,2,c(7,5))
#' ConvEncode(c(1,0,0,1,1), coder)
#' @author Martin Nocker
#' @export
ConvEncode <- function(message, conv.encoder, terminate = TRUE) {

  code <- c_convolutionEncode(message,
                              conv.encoder$N,
                              conv.encoder$M,
                              conv.encoder$next.state,
                              conv.encoder$output,
                              as.integer(conv.encoder$rsc),
                              conv.encoder$termination,
                              as.integer(terminate))
  return(code)
}

#' convolutional decoding of a code (soft decision)
#'
#' \code{ConvDecode} decodes a convolutional codeword
#' This decoder is a soft-input soft-output decoder
#' @param code the code to be decoded
#' @param conv.encoder convolutional encoder used for encoding [list]
#' @param terminate flag if the code was terminated, default: TRUE
#' @return the decoded message, list(softOutput, hardOutput)
#' @examples
#' coder <- GenerateConvEncoder(2,2,c(7,5))
#' coded <- ConvEncode(c(1,0,0,1,1), coder)
#' ConvDecode(coded, coder)
#' @author Martin Nocker
#' @export
ConvDecode <- function(code, conv.encoder, terminate = TRUE) {

  output <- c_convolutionDecode(code,
                                conv.encoder$N,
                                conv.encoder$M,
                                conv.encoder$prev.state,
                                conv.encoder$output)

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    M <- conv.encoder$M
    soft <- head(output$soft.output, length(output$soft.output) - M)
    hard <- head(output$hard.output, length(output$hard.output) - M)

    newlist <- list(soft.output = soft, hard.output = hard)
    return(newlist)
  }

  return(output)
}

#' convolutional decoding of a code (hard decision)
#'
#' \code{ConvDecodeHard} decodes a codeword that was encoded with the given
#' encoder. This decoder is a hard-decision decoder
#' @inheritParams ConvDecode
#' @return the hard-decoded message vector
#' @examples
#' coder <- GenerateConvEncoder(2,2,c(7,5))
#' coded <- ConvEncode(c(1,0,0,1,1), coder)
#' ConvDecodeHard(coded, coder)
#' @author Martin Nocker
#' @export
ConvDecodeHard <- function(code, conv.encoder, terminate = TRUE) {

  code.copy <- c(code)

  output <- c_convolutionDecode_hard(code.copy,
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
