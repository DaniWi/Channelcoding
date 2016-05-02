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

  generators <- generators[generators > 0]

  # encoder requires N generator polynoms
  if (length(generators) < N) {
    # stop execution if too few generators
    stop("Too few generator polynoms!")
  } else if (length(generators) > N) {
    # just use first N polynoms if too many are provided
    generators <- head(generators, N)
  }

  if (!IsOctal(generators)) {
    # only octal generators are accepted
    stop("At least one generator is not in octal form!")
  }

  max.generator.octal = DecimalToOctal(2^(M+1) - 1)

  if (any(generators > max.generator.octal)) {
    stop("At least one generator is greater than the maximum generator!")
    # generators = MaskGenerators(generators, max.generator.octal)
  }

  if (IsCatastrophicEncoder(generators)) {
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

  generators <- generators[generators > 0]

  # encoder requires N generator polynoms
  if (length(generators) < N) {
    # stop execution if too few generators
    stop("Too few generator polynoms!")
  } else if (length(generators) > N) {
    # just use first N polynoms if too many are provided
    generators <- head(generators, N)
  }

  if (!IsOctal(generators)) {
    # only octal generators are accepted
    stop("At least one generator is not in octal form!")
  }

  max.generator.octal = DecimalToOctal(2^(M+1) - 1)

  if (any(generators > max.generator.octal)) {
    stop("At least one generator is greater than the maximum generator!")
    # generators = MaskGenerators(generators, max.generator.octal)
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
ConvEncode <- function(message, conv.encoder = NULL, terminate = TRUE, punctuation.matrix = NULL, visualize = FALSE) {

  stopifnot(length(message) > 0)

  if (any((message != 1)[message != 0])) {
    stop("Nachricht darf nur 0er und 1er enthalten!")
  }

  if (is.null(conv.encoder)) {
    warning("Standard-Faltungskodierer wurde verwendet! N=2, M=2, Generatoren: (7,5)")
    conv.encoder <- GenerateConvEncoder(2,2,c(7,5))
  }

  if (!is.null(punctuation.matrix) && nrow(punctuation.matrix) != conv.encoder$N) {
    stop("Punktierungsmatrix hat falsche Anzahl an Zeilen! Matrix muss N Zeilen haben!")
  }

  code <- c_convolutionEncode(message,
                              conv.encoder$N,
                              conv.encoder$M,
                              conv.encoder$next.state,
                              conv.encoder$output,
                              as.integer(conv.encoder$rsc),
                              conv.encoder$termination,
                              as.integer(terminate))

  if (!is.null(punctuation.matrix)) {
    punctured.code <- PunctureCode(code, punctuation.matrix)

    return(list(original=code, punctured=punctured.code))
  }

  if (visualize) {
    rmarkdown::render(system.file("rmd", "ConvolutionEncode.Rmd", package = "channelcoding"),
                      encoding = "UTF-8",
                      params = list(conv.encoder = conv.encoder,
                                    message = message,
                                    terminate = terminate))

    rstudioapi::viewer(system.file("rmd", "ConvolutionEncode.pdf", package = "channelcoding"))
  }

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
ConvDecode <- function(code, conv.encoder = NULL, terminate = TRUE, punctuation.matrix = NULL, visualize = FALSE) {

  stopifnot(length(code) > 0)

  if (is.null(conv.encoder)) {
    warning("Standard-Faltungskodierer wurde verwendet! N=2, M=2, Generatoren: (7,5)")
    conv.encoder <- GenerateConvEncoder(2,2,c(7,5))
  }

  if(!is.null(punctuation.matrix)) {
    #insert missing bits from punctuation
    code <- InsertPunctuationBits(code, punctuation.matrix)
  }

  #Check length of input(with inserted bits) and permutation
  if ((length(code) %% conv.encoder$N) != 0) {
    if(!is.null(punctuation.matrix)) {
      stop("Fehler während der Punktierung!")
    }
    stop("Code hat die falsche Länge")
  }

  result <- c_convolutionDecode(code,
                                conv.encoder$N,
                                conv.encoder$M,
                                conv.encoder$prev.state,
                                conv.encoder$output,
                                as.integer(terminate))

  if (visualize) {
    rmarkdown::render(system.file("rmd", "ConvolutionDecode.Rmd", package = "channelcoding"),
                      encoding = "UTF-8",
                      params = list(conv.encoder = conv.encoder,
                                    code = code,
                                    decoded = result$hard.output,
                                    trellis = result$trellis,
                                    survivor.states = result$survivor.states,
                                    soft.flag = TRUE))

    rstudioapi::viewer(system.file("rmd", "ConvolutionDecode.pdf", package = "channelcoding"))
  }


  result <- result[1:2]

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    M <- conv.encoder$M
    soft <- head(result$soft.output, length(result$soft.output) - M)
    hard <- head(result$hard.output, length(result$hard.output) - M)

    newlist <- list(soft.output = soft, hard.output = hard)
    return(newlist)
  }

  return(result)
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
ConvDecodeHard <- function(code, conv.encoder, terminate = TRUE, punctuation.matrix = NULL) {

  stopifnot(length(code) > 0)

  if (is.null(conv.encoder)) {
    warning("Standard-Faltungskodierer wurde verwendet! N=2, M=2, Generatoren: (7,5)")
    conv.encoder <- GenerateConvEncoder(2,2,c(7,5))
  }

  code.copy <- c(code)

  if(!is.null(punctuation.matrix)) {
    #insert missing bits from punctuation
    code.copy <- InsertPunctuationBits(code.copy, punctuation.matrix)
  }

  #Check length of input(with inserted bits) and permutation
  if ((length(code.copy) %% conv.encoder$N) != 0) {
    if(!is.null(punctuation.matrix)) {
      stop("Fehler während der Punktierung!")
    }
    stop("Code hat die falsche Länge")
  }

  result <- c_convolutionDecode_hard(code.copy,
                                     conv.encoder$N,
                                     conv.encoder$M,
                                     conv.encoder$prev.state,
                                     conv.encoder$output,
                                     as.integer(terminate))

  rmarkdown::render(system.file("rmd", "ConvolutionDecode.Rmd", package = "channelcoding"),
                    encoding = "UTF-8",
                    params = list(conv.encoder = conv.encoder,
                                  code = code.copy,
                                  decoded = result$hard.output,
                                  trellis = result$trellis,
                                  survivor.states = result$survivor.states,
                                  soft.flag = FALSE))

  rstudioapi::viewer(system.file("rmd", "ConvolutionDecode.pdf", package = "channelcoding"))

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    return(head(result$hard.output, length(result$hard.output) - conv.encoder$M))
  }

  return(result$hard.output)
}

#' @export
ConvolutionSimulation <- function(coder = NULL,
                                  msg.length = 100,
                                  iterations.per.db = 100,
                                  min.db = 0.1,
                                  max.db = 3.0,
                                  db.interval = 0.1,
                                  punctuation.matrix = NULL)
{
  stopifnot(msg.length > 0, iterations.per.db > 0,
            min.db > 0, max.db > 0, max.db >= min.db, db.interval > 0)

  if (is.null(conv.encoder)) {
    warning("Standard-Faltungskodierer wurde verwendet! N=2, M=2, Generatoren: (7,5)")
    coder <- GenerateConvEncoder(2,2,c(7,5))
  }

  v.db <- seq(from = min.db, to = max.db, by = db.interval)
  v.ber <- numeric(0)

  total.errors <- 0

  for (db in v.db) {
    for (i in 1 : iterations.per.db) {
      # message erzeugen
      message <- sample(c(0,1), msg.length, replace = TRUE)

      # encode
      coded <- ConvEncode(message, coder, punctuation.matrix = punctuation.matrix)

      # wenn punktiert, muss codierte nachricht aus liste geholt werden
      if (!is.null(punctuation.matrix)) {
        coded <- coded$punctured
      }

      # noise hinzufügen
      noisy <- applyNoise(coded, db)

      # anzahl flipped bits (channel errors)
      # coded.hard <- ifelse(coded >= 0, 0, 1)
      # noisy.hard <- ifelse(noisy >= 0, 0, 1)
      # channel.errors <- sum(abs(coded.hard - noisy.hard))

      # decode
      decoded <- ConvDecode(noisy, coder, punctuation.matrix = punctuation.matrix)

      # vgl decoded & message
      decode.erros <- sum(abs(decoded$output.hard - message))

      total.errors <- total.errors + decode.erros
    }

    v.ber <- c(v.ber, total.errors / (msg.length * iterations.per.db))
    total.errors <- 0
  }


  df <- data.frame(db = v.db, ber = v.ber)

  return(df)
}
