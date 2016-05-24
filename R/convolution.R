# convolution.R
#
# interface for convolutional codes
# provided functions:
#  - generate non-recursive coder
#  - generate rsc coder
#  - encode
#  - decode (soft decision)
#  - decode (hard decision)
#  - simulation
#  - open pdf
#  - generate punctuation matrix

#' Generate convolutional encoder.
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
#' @param N Numer of output symbols per input symbol.
#' @param M Memory length of the encoder.
#' @param generators Vector of N octal generator polynoms
#'     (one for each output symbol).
#' @return A convolutional encoder represented as a list containing:
#'     N, M, vector of generator polynoms,
#'     3 matrices: nextState, previousState and output, rsc flag,
#'     termination vector
#' @examples
#' # standard convolutional encoder with code-rate = 0.5
#' ConvGenerateEncoder(2,2,c(7,5))
#' @author Martin Nocker
#' @export
ConvGenerateEncoder <- function(N, M, generators) {

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

  if (IsCatastrophicEncoder(generators, M+1)) {
    warning("The generated encoder is a catastrophic encoder!")
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

#' Generate rsc encoder.
#'
#' Generates a recursive systematic convolutional (rsc) encoder.
#' @details N is an integer and gives the number of output bits per input bit.
#'     N has to be at least two. M is an integer and gives the memory length
#'     of the encoder (number of shift register elements in the circuit). M
#'     has to be at least one.
#'     The generator polynoms define how the output bits are computed for each
#'     of the N output signals. The polynoms are octal numbers. See details of
#'     \code{\link{ConvGenerateEncoder}} for an example.
#'     An rsc encoder has exactly one fixed systematic output signal.
#'     The generator polynom for the systematic output doesn't have to be
#'     passed as an argument. So the generators argument contains all polynoms
#'     for non-systematic outputs and at the last position the recursion
#'     polynom. The MSB of the recursion polynom handles the input signal,
#'     the other bits handle the memory outputs. The MSB of the output polynoms
#'     handle the recursion output(!), not the original input signal. The other
#'     bits also handle the memory outputs. See 'The art of error correcting
#'     coding' (p.92f) for a detailed definition and an example.
#' @param N Numer of output symbols per input symbol.
#' @param M Memory length of the encoder.
#' @param generators Vector of octal generator polynoms
#'     (one for each non-systematic output symbol and one for the recursion).
#' @return A convolutional encoder represented as a list containing:
#'     N, M, vector of generator polynoms,
#'     3 matrices: nextState, previousState and output, rsc flag,
#'     termination vector
#' @examples
#' # standard rsc encoder with code-rate = 0.5
#' ConvGenerateRscEncoder(2,2,c(5,7))
#' @author Martin Nocker
#' @export
ConvGenerateRscEncoder <- function(N, M, generators) {

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

#' Convolutional encoding of a message.
#'
#' Produces a convolutional code of a message.
#' @param message The message to be encoded.
#' @param conv.encoder Convolutional encoder used for encoding.
#' @param terminate Flag if the code should be terminated.
#' @param punctuation.matrix If not null the encoded message is punctured with the punctuation matrix.
#' @param visualize If TRUE a beamer PDF file is generated showing the encode process.
#' @return The encoded message with signal values {-1,+1} which map to {1,0} respectively.
#' @examples
#' coder <- ConvGenerateEncoder(2,2,c(7,5))
#' plain <- c(1,0,0,1,1)
#'
#' # standard encoding with termination
#' ConvEncode(plain, coder)
#'
#' # without termination
#' ConvEncode(plain, coder, terminate = FALSE)
#'
#' # with punctuation and visualization
#' # punctuation <- ConvGetPunctuationMatrix(c(1,1,0,1,1,0), coder)
#' # ConvEncode(plain, coder, punctuation.matrix = punctuation, visualize = TRUE)
#'
#' # use default values
#' ConvEncode(plain)
#' @author Martin Nocker
#' @export
ConvEncode <- function(message,
                       conv.encoder = NULL,
                       terminate = TRUE,
                       punctuation.matrix = NULL,
                       visualize = FALSE)
{
  stopifnot(length(message) > 0)

  if (any((message != 1)[message != 0])) {
    stop("Message should only contain 0 and 1!")
  }

  if (is.null(conv.encoder)) {
    warning("Standard Convolutional Coder was used! N=2, M=2, Generators: (7,5)")
    conv.encoder <- ConvGenerateEncoder(2,2,c(7,5))
  }

  CheckCoder(conv.encoder)

  if (!is.null(punctuation.matrix) && nrow(punctuation.matrix) != conv.encoder$N) {
    stop("Punctuation matrix has the wrong amount of rows! Matirx has to have N rows!")
  }

  code <- c_convolutionEncode(message,
                              conv.encoder$N,
                              conv.encoder$M,
                              conv.encoder$next.state,
                              conv.encoder$output,
                              as.integer(conv.encoder$rsc),
                              conv.encoder$termination,
                              as.integer(terminate))

  if (visualize) {
    if ((length(message) + as.integer(terminate)*conv.encoder$M) > 14) {
      warning("Message is too long for a proper visualization! No PDF was created! Maximum length is 14 (including termination)")
    } else if (is.null(punctuation.matrix)) {
      if (conv.encoder$M > 3) {
        warning("Coder has more than 8 states (M > 3). A PDF was generated but without a coder visualization!")
      }
      rmarkdown::render(system.file("rmd", "ConvolutionEncode.Rmd", package = "channelcoding"),
                        output_dir = system.file("pdf", package = "channelcoding"),
                        encoding = "UTF-8",
                        params = list(conv.encoder = conv.encoder,
                                      message = message,
                                      code = code,
                                      terminate = terminate))

      rstudioapi::viewer(system.file("pdf", "ConvolutionEncode.pdf", package = "channelcoding"))
    } else {
      if (conv.encoder$M > 3) {
        warning("Coder has more than 8 states (M > 3). A PDF was generated but without a coder visualization!")
      }
      rmarkdown::render(system.file("rmd", "ConvolutionEncodePunctured.Rmd", package = "channelcoding"),
                        output_dir = system.file("pdf", package = "channelcoding"),
                        encoding = "UTF-8",
                        params = list(conv.encoder = conv.encoder,
                                      message = message,
                                      code = code,
                                      terminate = terminate,
                                      punctuation = punctuation.matrix))

      rstudioapi::viewer(system.file("pdf", "ConvolutionEncodePunctured.pdf", package = "channelcoding"))
    }
  }

  if (!is.null(punctuation.matrix)) {
    punctured.code <- PunctureCode(code, punctuation.matrix)

    return(list(original=code, punctured=punctured.code))
  }

  return(code)
}

#' Convolutional decoding of a code (soft decision).
#'
#' Decodes a convolutional codeword using soft decision decoding.
#' @param code The code to be decoded.
#' @param conv.encoder Convolutional encoder used for encoding.
#' @param terminate flag If the code is terminated.
#' @param punctuation.matrix If not null the code is depunctured prior to the decode algorithm.
#' @param visualize If TRUE a beamer PDF file is generated showing the decode process.
#' @return The decoded message, list(softOutput, hardOutput)
#' @examples
#' coder <- ConvGenerateEncoder(2,2,c(7,5))
#' plain <- c(1,0,0,1,1)
#'
#' # standard encoding and decoding with termination
#' coded <- ConvEncode(plain, coder)
#' ConvDecodeSoft(coded, coder)
#'
#' # without termination
#' coded <- ConvEncode(plain, coder, terminate = FALSE)
#' ConvDecodeSoft(coded, coder, terminate = FALSE)
#'
#' # with punctuation and visualization
#' # punctuation <- ConvGetPunctuationMatrix(c(1,1,0,1,1,0), coder)
#' # coded <- ConvEncode(plain, coder, punctuation.matrix = punctuation)
#' # ConvDecodeSoft(coded, coder, punctuation.matrix = punctuation, visualize = TRUE)
#'
#' # with message distortion
#' coded <- ConvEncode(plain, coder)
#' noisy <- ApplyNoise(coded, SNR.db = 3)
#' ConvDecodeSoft(noisy, coder)
#' @author Martin Nocker
#' @export
ConvDecodeSoft <- function(code,
                           conv.encoder = NULL,
                           terminate = TRUE,
                           punctuation.matrix = NULL,
                           visualize = FALSE)
{
  if (!is.numeric(code) || !is.vector(code)) {
    stop("Code was not a numeric vector! Check the code parameter!")
  }

  stopifnot(length(code) > 0)

  if (is.null(conv.encoder)) {
    warning("Standard Convolutional Coder was used! N=2, M=2, Generators: (7,5)")
    conv.encoder <- ConvGenerateEncoder(2,2,c(7,5))
  }

  CheckCoder(conv.encoder)

  punct.code <- 0  # needed for visualization in case code in punctured
  if(!is.null(punctuation.matrix)) {
    punct.code <- c(code)
    #insert missing bits from punctuation
    code <- InsertPunctuationBits(code, punctuation.matrix)
  }

  #Check length of input(with inserted bits) and permutation
  if ((length(code) %% conv.encoder$N) != 0) {
    if(!is.null(punctuation.matrix)) {
      stop("Mistake during depunctuation!")
    }
    stop("Code has the wrong length! (Code is maybe punctured)")
  }

  result <- c_convolutionDecode(code,
                                conv.encoder$N,
                                conv.encoder$M,
                                conv.encoder$prev.state,
                                conv.encoder$output,
                                as.integer(terminate))

  if (visualize) {
    if (conv.encoder$M > 3) {
      warning("Coder has more than 8 states (M > 3) and thus can't be visualized! No PDF was created!")
    } else if (length(result$output.hard) > 14) {
      warning("Code is too long for a proper visualization! No PDF was created! Maximum length is 14 (including termination)")
    } else if (is.null(punctuation.matrix)) {
      rmarkdown::render(system.file("rmd", "ConvolutionDecode.Rmd", package = "channelcoding"),
                        output_dir = system.file("pdf", package = "channelcoding"),
                        encoding = "UTF-8",
                        params = list(conv.encoder = conv.encoder,
                                      code = code,
                                      decoded = result$output.hard,
                                      trellis = result$trellis,
                                      survivor.states = result$survivor.states,
                                      soft.flag = TRUE))

      rstudioapi::viewer(system.file("pdf", "ConvolutionDecode.pdf", package = "channelcoding"))
    } else {
      rmarkdown::render(system.file("rmd", "ConvolutionDecodePunctured.Rmd", package = "channelcoding"),
                        output_dir = system.file("pdf", package = "channelcoding"),
                        encoding = "UTF-8",
                        params = list(conv.encoder = conv.encoder,
                                      code = code,
                                      decoded = result$output.hard,
                                      trellis = result$trellis,
                                      survivor.states = result$survivor.states,
                                      punctuation = punctuation.matrix,
                                      punctured.code = punct.code,
                                      soft.flag = TRUE))
      rstudioapi::viewer(system.file("pdf", "ConvolutionDecodePunctured.pdf", package = "channelcoding"))
    }
  }


  result <- result[1:2]

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    M <- conv.encoder$M
    soft <- head(result$output.soft, length(result$output.soft) - M)
    hard <- head(result$output.hard, length(result$output.hard) - M)

    newlist <- list(output.soft = soft, output.hard = hard)
    return(newlist)
  }

  return(result)
}

#' Convolutional decoding of a code (hard decision).
#'
#' Decodes a convolutional codeword using hard decision decoding.
#' @inheritParams ConvDecodeSoft
#' @return The hard-decoded message vector.
#' @examples
#' coder <- ConvGenerateEncoder(2,2,c(7,5))
#' plain <- c(1,0,0,1,1)
#'
#' # standard encoding and decoding with termination
#' coded <- ConvEncode(plain, coder)
#' ConvDecodeHard(coded, coder)
#'
#' # without termination
#' coded <- ConvEncode(plain, coder, terminate = FALSE)
#' ConvDecodeHard(coded, coder, terminate = FALSE)
#'
#' # with punctuation and visualization
#' # punctuation <- ConvGetPunctuationMatrix(c(1,1,0,1,1,0), coder)
#' # coded <- ConvEncode(plain, coder, punctuation.matrix = punctuation)
#' # ConvDecodeHard(coded, coder, punctuation.matrix = punctuation, visualize = TRUE)
#'
#' # with message distortion
#' coded <- ConvEncode(plain, coder)
#' noisy <- ApplyNoise(coded, SNR.db = 3)
#' ConvDecodeHard(noisy, coder)
#' @author Martin Nocker
#' @export
ConvDecodeHard <- function(code,
                           conv.encoder = NULL,
                           terminate = TRUE,
                           punctuation.matrix = NULL,
                           visualize = FALSE)
{
  if (!is.numeric(code) || !is.vector(code)) {
    stop("Code was not a numeric vector! Check the code parameter!")
  }

  stopifnot(length(code) > 0)

  if (is.null(conv.encoder)) {
    warning("Standard Convolutional Coder was used! N=2, M=2, Generators: (7,5)")
    conv.encoder <- ConvGenerateEncoder(2,2,c(7,5))
  }

  CheckCoder(conv.encoder)

  code.copy <- c(code)

  punct.code <- 0 # needed for visualization in case code in punctured
  if(!is.null(punctuation.matrix)) {
    punct.code <- c(code.copy)
    #insert missing bits from punctuation
    code.copy <- InsertPunctuationBits(code.copy, punctuation.matrix)
  }

  #Check length of input(with inserted bits) and permutation
  if ((length(code.copy) %% conv.encoder$N) != 0) {
    if(!is.null(punctuation.matrix)) {
      stop("Mistake during depunctuation!")
    }
    stop("Code has the wrong length! (Code is maybe punctured)")
  }

  result <- c_convolutionDecode_hard(code.copy,
                                     conv.encoder$N,
                                     conv.encoder$M,
                                     conv.encoder$prev.state,
                                     conv.encoder$output,
                                     as.integer(terminate))

  if (visualize) {
    if (conv.encoder$M > 3) {
      warning("Coder has more than 8 states (M > 3) and thus can't be visualized! No PDF was created!")
    } else if (length(result$output.hard) > 14) {
      warning("Code is too long for a proper visualization! No PDF was created! Maximum length is 14 (including termination)")
    } else if (is.null(punctuation.matrix)) {
      rmarkdown::render(system.file("rmd", "ConvolutionDecode.Rmd", package = "channelcoding"),
                        output_dir = system.file("pdf", package = "channelcoding"),
                        encoding = "UTF-8",
                        params = list(conv.encoder = conv.encoder,
                                      code = code,
                                      decoded = result$output.hard,
                                      trellis = result$trellis,
                                      survivor.states = result$survivor.states,
                                      soft.flag = FALSE))

      rstudioapi::viewer(system.file("pdf", "ConvolutionDecode.pdf", package = "channelcoding"))
    } else {
      rmarkdown::render(system.file("rmd", "ConvolutionDecodePunctured.Rmd", package = "channelcoding"),
                        output_dir = system.file("pdf", package = "channelcoding"),
                        encoding = "UTF-8",
                        params = list(conv.encoder = conv.encoder,
                                      code = code.copy,
                                      decoded = result$output.hard,
                                      trellis = result$trellis,
                                      survivor.states = result$survivor.states,
                                      punctuation = punctuation.matrix,
                                      punctured.code = punct.code,
                                      soft.flag = FALSE))
      rstudioapi::viewer(system.file("pdf", "ConvolutionDecodePunctured.pdf", package = "channelcoding"))
    }
  }

  # if terminated, termination bits are thrown away
  if (terminate == TRUE) {
    return(head(result$output.hard, length(result$output.hard) - conv.encoder$M))
  }

  return(result$output.hard)
}

#' Convolutional Simulation.
#'
#' Simulation of a convolutional encode and decode process over a noisy channel
#' at several signal-noise-ratios (SNR).
#'
#' @param conv.coder Convolutional coder used for the simulation. Can be created via
#'     \code{\link{ConvGenerateEncoder}} or \code{\link{ConvGenerateRscEncoder}}.
#' @param msg.length Message length of the randomly created messages to be encoded.
#' @param min.db Minimum SNR to be tested.
#' @param max.db Maximum SNR to be tested.
#' @param db.interval Step between two SNRs tested.
#' @param iterations.per.db Number of encode and decode processes per SNR.
#' @param punctuation.matrix If not null the process involves the punctuation. Can
#'     be created via \code{\link{ConvGetPunctuationMatrix}}.
#' @param visualize If true a PDF report is generated.
#' @return Dataframe containing the bit-error-rates for each SNR tested.
#' @examples
#' # use all default parameters
#' ConvSimulation()
#'
#' # without punctuation
#' coder <- ConvGenerateEncoder(2,2,c(7,5))
#' ConvSimulation(coder, 10, 0.01, 1, 0.05, 50, NULL, FALSE)
#' @author Martin Nocker
#' @export
ConvSimulation <- function(conv.coder = NULL,
                           msg.length = 100,
                           min.db = 0.1,
                           max.db = 2.0,
                           db.interval = 0.1,
                           iterations.per.db = 100,
                           punctuation.matrix = NULL,
                           visualize = FALSE)
{
  stopifnot(msg.length > 0, iterations.per.db > 0,
            min.db > 0, max.db > 0, max.db >= min.db, db.interval > 0)

  if (is.null(conv.coder)) {
    warning("Standard Convolutional Coder was used! N=2, M=2, Generators: (7,5)")
    conv.coder <- ConvGenerateEncoder(2,2,c(7,5))
  }

  CheckCoder(conv.coder)

  v.db <- seq(from = min.db, to = max.db, by = db.interval)
  v.ber <- numeric(0)

  total.errors <- 0

  total.iterations <- length(v.db) * iterations.per.db

  print("Convolution Simulation")
  progress.bar <- txtProgressBar(min = 0, max = total.iterations, style = 3)

  progress.counter <- 0

  for (db in v.db) {
    for (i in 1 : iterations.per.db) {
      # message erzeugen
      message <- sample(c(0,1), msg.length, replace = TRUE)

      # encode
      coded <- ConvEncode(message, conv.coder, punctuation.matrix = punctuation.matrix)

      # wenn punktiert, muss codierte nachricht aus liste geholt werden
      if (!is.null(punctuation.matrix)) {
        coded <- coded$punctured
      }

      # noise hinzufügen
      noisy <- ApplyNoise(coded, db)

      # anzahl flipped bits (channel errors)
      # coded.hard <- ifelse(coded >= 0, 0, 1)
      # noisy.hard <- ifelse(noisy >= 0, 0, 1)
      # channel.errors <- sum(abs(coded.hard - noisy.hard))

      # decode
      decoded <- ConvDecodeSoft(noisy, conv.coder, punctuation.matrix = punctuation.matrix)

      # vgl decoded & message
      decode.errors <- sum(abs(decoded$output.hard - message))

      total.errors <- total.errors + decode.errors

      progress.counter <- progress.counter + 1

      setTxtProgressBar(progress.bar, progress.counter)
    }

    v.ber <- c(v.ber, total.errors / (msg.length * iterations.per.db))
    total.errors <- 0
  }

  close(progress.bar)

  df <- data.frame(db = v.db, ber = v.ber)

  if (visualize) {
    rmarkdown::render(system.file("rmd", "Simulation.Rmd", package = "channelcoding"),
                      output_dir = system.file("pdf", package = "channelcoding"),
                      output_file = "ConvolutionSimulation.pdf",
                      encoding = "UTF-8",
                      params = list(turbo = FALSE,
                                    message.length = msg.length,
                                    iterations.per.db = iterations.per.db,
                                    min.db = min.db,
                                    max.db = max.db,
                                    db.interval = db.interval,
                                    punctuation = punctuation.matrix,
                                    encoder = conv.coder,
                                    dataframe = df))

    rstudioapi::viewer(system.file("pdf", "ConvolutionSimulation.pdf", package = "channelcoding"))
  }

  return(df)
}

#' Open visualization PDF.
#'
#' With this function it is easy to reopen the PDF files created with
#'     \code{\link{ConvEncode}}, \code{\link{ConvDecodeSoft}},
#'     \code{\link{ConvDecodeHard}} and \code{\link{ConvSimulation}}.
#'     The files are stored in the program files of R.
#'     If the corresponding file does not exist yet an error message is printed.
#' @param encode Flag to open encode PDFs (if TRUE) or decode PDFs (if FALSE).
#' @param punctured Flag to open encode or decode PDFs with punctuation.
#' @param simulation Flag to open simulation PDFs. This flag has highest precedence
#'     meaning that if the simulation flag is TRUE the function will look for the
#'     simulation PDF. Only if FALSE the others are evaluated.
#' @examples
#' # open encode without punctuation PDF
#' ConvOpenPDF(encode = TRUE, punctured = FALSE)
#'
#' # open decode with punctuation PDF
#' ConvOpenPDF(encode = FALSE, punctured = TRUE)
#'
#' # open convolution simulation PDF
#' ConvOpenPDF(simulation = TRUE)
#' @author Martin Nocker
#' @export
ConvOpenPDF <- function(encode = TRUE, punctured = FALSE, simulation = FALSE) {
  if (simulation) {
    path <- system.file("pdf", "ConvolutionSimulation.pdf", package = "channelcoding")
  } else {
    if (encode) {
      if (punctured) {
        path <- system.file("pdf", "ConvolutionEncodePunctured.pdf", package = "channelcoding")
      } else {
        path <- system.file("pdf", "ConvolutionEncode.pdf", package = "channelcoding")
      }
    } else {
      if (punctured) {
        path <- system.file("pdf", "ConvolutionDecodePunctured.pdf", package = "channelcoding")
      } else {
        path <- system.file("pdf", "ConvolutionDecode.pdf", package = "channelcoding")
      }
    }
  }
  if (path != "") {
    rstudioapi::viewer(path)
  } else {
    warning("File does not exists!")
  }
}

#' Get punctuation matrix from vector.
#'
#' Creates a punctuation matrix from the passed punctuation vector and the
#' passed coder.
#' Important: punctuation vector has to contain k*N elements, where N is the
#' amount of coder outputs and k a natural number \eqn{\ge 1}.
#' @param punctuation.vector Vector containing the punctuation information which will
#'     be transformed to a punctuation matrix.
#' @param conv.coder Convolutional coder which is used for the matrix dimensions.
#' @return Punctuation matrix suitable for \code{\link{ConvEncode}},
#'     \code{\link{ConvDecodeSoft}}, \code{\link{ConvDecodeHard}} and
#'     \code{\link{ConvSimulation}}.
#' @examples
#' # coder with 2 outputs, vector has to have k*2 elements
#' coder <- ConvGenerateEncoder(2,2,c(7,5))
#' punctuation.vector <- c(1,1,0,1,1,0)
#' punctuation.matrix <- ConvGetPunctuationMatrix(punctuation.vector, coder)
#' @author Martin Nocker
#' @export
ConvGetPunctuationMatrix <- function(punctuation.vector, conv.coder) {

  if (any((punctuation.vector != 1)[punctuation.vector != 0])) {
    # punctuation.vector has elements with value different from 0 or 1 which is not allowed
    stop("Invalid punctuation vector! Only values 0/1 are allowed!")
  }

  if (is.null(conv.coder$N)) {
    stop("Encoder has to specify list element N!")
  }

  if (length(punctuation.vector) %% conv.coder$N != 0) {
    stop("Wrong length of punctuation vector! Must be a multiple of N (amount of exits)!")
  }

  mat <- matrix(punctuation.vector, nrow = conv.coder$N)

  if (any(colSums(mat) == 0)) {
    stop("Punctuation matrix should not have a 0 column!")
  }

  return(mat)
}
