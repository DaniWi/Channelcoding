# turbo_code_rcpp.R
#
# interface for turbo codes
# provided functions:
#  - TurboEncode
#  - TurboDecode
#  - TurboGetpermutation
#  - TurboGetPunctuationMatrix
#  - TurboSimutlation
#  - TurboOpenPDF


#' Encode a message with the turbo-code.
#'
#' This functions takes a message and encode it with the turbo-code procedure. The message
#' will be encoded in 2 identical systematic encoder, but before the message will be encoded in the encoder 2
#' the message will be permutated with an interleaver. To increase the minimal distance between 2 code words
#' the interleaver is needed. The systematic encoder must have 2 exists and it must be
#' an systematic encoder which means that the input is redirected directly to the output.
#' Because of the systematic encoders the original message must be transmitted only once, so
#' this procedure reduces the code rate of the whole encoder. After the 2 encoders the
#' resulting encoded messages are interleaved into one stream and transmitted to the exit.
#' Only one exit of the systematic encoder will be used, so the index of the output can be
#' determined. The punctuation matrix removes all bits where the matrix is zero, so the
#' code rate is increased.
#'
#' @author Witsch Daniel
#'
#' @param message Message which will be encoded.
#' @param permutation.vector Permutation vector which will be created with \code{\link{TurboGetPermutation}}.
#' @param coder.info Coder which will be created with \code{\link{ConvGenerateEncoder}} or \code{\link{ConvGenerateRscEncoder}}.
#' @param parity.index Index to decide which exit of the coder will be used for the encoding (>1).
#'        Default value is always the last exit of the coder.
#' @param punctuation.matrix Punctuation matrix to puncture the output, will be created with \code{\link{TurboGetPunctuationMatrix}}.
#' @param visualize Flag to decide whether to create a visualization pdf or not.
#'
#' @return Encoded message, will be a list if the output was punctured (list with original and punctured).
#'
#' @examples
#' input <- c(1,0,1,1,0)
#'
#' #default coder and permutation vector
#' message.encoded <- TurboEncode(input)
#' print(message.encoded)
#'
#' #custom coder and permutation vector
#' coder <- ConvGenerateRscEncoder(2, 2, c(5, 7))
#' perm <- TurboGetPermutation(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#' print(message.encoded)
#'
#' #encoding with punctuation matrix
#' punct <- TurboGetPunctuationMatrix(c(1,1,0,1,0,1))
#' message.encoded.punct <- TurboEncode(input, perm, coder, punctuation.matrix = punct)
#' print(message.encoded.punct$original)
#' print(message.encoded.punct$punctured)
#'
#' @export
TurboEncode <-
  function(message,
           permutation.vector = NULL,
           coder.info = NULL,
           parity.index = coder.info$N,
           punctuation.matrix = NULL,
           visualize = FALSE) {

    if (is.null(coder.info)) {
      warning("Standard-Coder was used! RSC, N=2, M=2, Generators: (5,7)")
      coder.info <- ConvGenerateRscEncoder(2,2,c(5,7))
    }
    if (is.null(permutation.vector)) {
      warning("Standard-Permutationvector was used! (PRIMITIVE)")
      permutation.vector <- TurboGetPermutation(length(message), coder.info, "PRIMITIVE", list(root=0))
    }

    #Checks all parameters
    stopifnot(length(message) > 0)

    CheckCoder(coder.info)

    if (any((message != 1)[message != 0])) {
      stop("Message should only contain 0 and 1!")
    }
    if (coder.info$N < 2) {
      stop("Coder must have minimum 2 exists!")
    }
    if (coder.info$rsc != TRUE && coder.info$generators[1] != (2 ^ coder.info$M)) {
      stop("Coder must be an systematic coder!
            (exit 1 must have polynom 2^M or it must be a rescursiv coder)")
    }
    if (length(permutation.vector) != (length(message) + coder.info$M)) {
      stop("Permutation vector has the wrong length!")
    }
    if (any(is.na(match(0:(length(permutation.vector) - 1), permutation.vector)))) {
      stop("Permutation vector has the wrong entries!")
    }
    if (!is.null(punctuation.matrix) && nrow(punctuation.matrix) != 3) {
      stop("Punctuation matrix has the wrong amount of rows,
            with turbo codes it has to be 3!")
    }

    #when parity.index is too high or to low
    if (parity.index > coder.info$N || parity.index < 0) {
      parity.index = coder.info$N
    }

    #first coder with termination
    parity.1 <- ConvEncode(message, coder.info, TRUE)
    #extract original message with termination bits from parity
    temp.index <- c(rep(FALSE, 0), TRUE, rep(FALSE, coder.info$N - 1))
    code.orig <- parity.1[temp.index]

    #permutate original message with termination (mapping important!!)
    code.orig.01 <- as.numeric(code.orig < 0)
    code.perm <- as.numeric(code.orig[permutation.vector + 1] < 0)

    #second coder without termination
    parity.2 <- ConvEncode(code.perm, coder.info, FALSE)


    #extract parity bits from the output messages
    temp.index <- c(rep(FALSE, parity.index - 1), TRUE, rep(FALSE, coder.info$N - parity.index))
    parity.1 <- parity.1[temp.index]
    parity.2 <- parity.2[temp.index]

    code.result <- Interleave(code.orig, parity.1, parity.2)

    if (!is.null(punctuation.matrix)) {
      #puncture the output
      code.punct <- PunctureCode(code.result, punctuation.matrix)
    }

    if (visualize) {
      if (length(message) > 10) {
        warning("Maybe not the whole informations could be visualized!")
      }
      if (!is.null(punctuation.matrix)) {
        rmarkdown::render(
          system.file("rmd", "TurboEncodePunctured.Rmd", package = "channelcoding"),
          output_dir = system.file("pdf", package = "channelcoding"),
          params = list(
            orig = code.orig.01,
            interl = code.perm,
            parity1 = parity.1,
            parity2 = parity.2,
            multipl = code.result,
            result = code.punct,
            permutation = permutation.vector,
            punctuation = punctuation.matrix,
            encoder = coder.info),
          encoding = "UTF-8")
        rstudioapi::viewer(system.file("pdf", "TurboEncodePunctured.pdf", package = "channelcoding"))
      } else {
        rmarkdown::render(
          system.file("rmd", "TurboEncode.Rmd", package = "channelcoding"),
          output_dir = system.file("pdf", package = "channelcoding"),
          params = list(
            orig = code.orig.01,
            interl = code.perm,
            parity1 = parity.1,
            parity2 = parity.2,
            result = code.result,
            permutation = permutation.vector,
            encoder = coder.info),
          encoding = "UTF-8")
        rstudioapi::viewer(system.file("pdf", "TurboEncode.pdf", package = "channelcoding"))
      }
    }

    if (!is.null(punctuation.matrix)) {
      return(list(original = code.result, punctured = code.punct))
    } else {
      return(code.result)
    }

  }


#' Decode a code which is encoded with the turbo-code.
#'
#' The code will be decoded with the turbo-code procedure. The given coder must be an
#' systematic coder and should have at least 2 exits. To decide which exit is used, the user
#' can decide this with the parameter index. The principle of turbo decoding is that
#' the decoding is a iterative method, so code will be passed through the decoder multiple times.
#' With the parameter iteration you can decide how often the decoders are passed through.
#' After each iteration the result should be better and the result of the last iteration will
#' be used in the next iteration as input values. The exact signal level of the input will be used
#' as input in the first iteration. When punctuation was used during encoding it is important
#' that the same punctuation is given in the decoding because all bits which are deleted can now
#' be inserted again.
#'
#' @author Witsch Daniel
#'
#' @param code Code which will be decoded to the original message.
#' @param permutation.vector Permutation vector which will be created with \code{\link{TurboGetPermutation}}.
#' @param iterations Amount of decoding iterations.
#' @param coder.info Coder which will be created with \code{\link{ConvGenerateEncoder}} or \code{\link{ConvGenerateRscEncoder}}.
#' @param parity.index Index to decide which exit of the coder will be used for the encoding (>1).
#'        Default value is always the last exit of the coder.
#' @param punctuation.matrix Punctuation matrix to puncture the output, will be created with \code{\link{TurboGetPunctuationMatrix}}.
#' @param visualize Flag to decide whether to create a visualization pdf or not.
#'
#' @return Decoded message, will be a list (soft and hard values).
#'
#' @examples
#' input <- c(1,0,1,1,0)
#'
#' #default coder and permutation vector
#' message.encoded <- TurboEncode(input)
#' result <- TurboDecode(message.encoded, iterations = 5)
#' print(result)
#'
#' #custom coder and permutation vector
#' coder <- ConvGenerateRscEncoder(2, 2, c(5, 7))
#' perm <- TurboGetPermutation(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#' result <- TurboDecode(message.encoded, perm, 5, coder)
#' print(result)
#'
#' #decoding with punctuation matrix
#' punct <- TurboGetPunctuationMatrix(c(1,1,0,1,0,1))
#' message.encoded.punct <- TurboEncode(input, perm, coder, punctuation.matrix = punct)
#' result <- TurboDecode(message.encoded.punct$punctured, perm, 5, coder, punctuation.matrix = punct)
#' print(result)
#'
#' @export
TurboDecode <-
  function(code,
           permutation.vector = NULL,
           iterations = 1,
           coder.info = NULL,
           parity.index = coder.info$N,
           punctuation.matrix = NULL,
           visualize = FALSE) {

    if (!is.numeric(code) || !is.vector(code)) {
      stop("Code was not a numeric vector! Check the code parameter!")
    }

    if (is.null(coder.info)) {
      warning("Standard-Coder was used! RSC, N=2, M=2, Generators: (5,7)")
      coder.info <- ConvGenerateRscEncoder(2,2,c(5,7))
    }

    #Checks all parameters
    stopifnot(length(code) > 0, iterations > 0)

    CheckCoder(coder.info)

    if (coder.info$N < 2) {
      stop("Coder must have minimum 2 exists!")
    }
    if (coder.info$rsc != TRUE && coder.info$generators[1] != (2 ^ coder.info$M)) {
      stop("Coder must be an systematic coder!
            (exit 1 must have polynom 2^M or it must be a rescursiv coder)")
    }
    if (!is.null(permutation.vector) && any(is.na(match(0:(length(permutation.vector) - 1), permutation.vector)))) {
      stop("Permutation vector has the wrong entries!")
    }
    if (!is.null(punctuation.matrix) && nrow(punctuation.matrix) != 3) {
      stop("Punctuation matrix has the wrong amount of rows,
            with turbo codes it has to be 3!")
    }

    #when parity.index is too high or to low
    if (parity.index > coder.info$N || parity.index < 0) {
      parity.index = coder.info$N
    }

    #Punctuation
    if (!is.null(punctuation.matrix)) {
      #insert missing bits from punctuation
      code.with.punct <- InsertPunctuationBits(code, punctuation.matrix)
    } else {
      code.with.punct <- code
    }

    #create the default permutation.vector
    if (is.null(permutation.vector)) {
      warning("Standard-Permutationvector was used! (PRIMITIVE)")
      permutation.vector <- TurboGetPermutation(length(code.with.punct) / 3  - coder.info$M, coder.info, "PRIMITIVE", list(root=0))
    }

    #Check length of input(with inserted bits) and permutation.vector
    if ((length(code.with.punct) %% 3) != 0) {
      if (!is.null(punctuation.matrix)) {
        stop("Mistake during depunctation!")
      }
      stop("Code has the wrong length! (Code is maybe punctuated)")
    }
    if (length(permutation.vector) != (length(code.with.punct) / 3)) {
      stop("Permutation vector has the wrong length!")
    }

    #parse original parts of the code from input code
    code.length <- length(code.with.punct) / 3
    code.orig <- Deinterleave(code.with.punct, 1)
    parity.1 <- Deinterleave(code.with.punct, 2)
    parity.2 <- Deinterleave(code.with.punct, 3)

    #decode message (c-function)
    decoded <-
      c_turbo_decode(
        code.orig,
        parity.1,
        parity.2,
        permutation.vector,
        iterations,
        coder.info$N,
        coder.info$M,
        coder.info$prev.state,
        coder.info$output,
        parity.index
      )

    #delete termination bits from result
    output.soft <- head(decoded$soft.output, code.length - coder.info$M)
    output.hard <- head(decoded$hard.output, code.length - coder.info$M)
    message.decoded <- list(output.soft = output.soft, output.hard = output.hard)

    if (visualize) {
      if (length(output.hard) > 10 || iterations > 10) {
        warning("Maybe not the whole informations could be visualized!")
      }
      if (!is.null(punctuation.matrix)) {
        rmarkdown::render(
          system.file("rmd", "TurboDecodePunctured.Rmd", package = "channelcoding"),
          output_dir = system.file("pdf", package = "channelcoding"),
          params = list(
            code = code,
            punctured = code.with.punct,
            orig = code.orig,
            parity1 = parity.1,
            parity2 = parity.2,
            origI = decoded$disp.info$origI,
            decode1 = decoded$disp.info$decode1,
            decode1I = decoded$disp.info$decode1I,
            decode2Back = head(decoded$disp.info$decode2Back, iterations - 1),
            decode2IBack = head(decoded$disp.info$decode2IBack, iterations - 1),
            decode2 = tail(decoded$disp.info$decode2Back, 1)[[1]],
            decode2I = tail(decoded$disp.info$decode2IBack, 1)[[1]],
            tempResultSoft = decoded$disp.info$tempResultSoft,
            tempResultHard = decoded$disp.info$tempResultHard,
            result = message.decoded,
            iterations = iterations,
            permutation = permutation.vector,
            punctuation = punctuation.matrix,
            encoder = coder.info),
          encoding = "UTF-8")
        rstudioapi::viewer(system.file("pdf", "TurboDecodePunctured.pdf", package = "channelcoding"))
      } else {
        rmarkdown::render(
          system.file("rmd", "TurboDecode.Rmd", package = "channelcoding"),
          output_dir = system.file("pdf", package = "channelcoding"),
          params = list(
            code = code,
            orig = code.orig,
            parity1 = parity.1,
            parity2 = parity.2,
            origI = decoded$disp.info$origI,
            decode1 = decoded$disp.info$decode1,
            decode1I = decoded$disp.info$decode1I,
            decode2Back = head(decoded$disp.info$decode2Back, iterations - 1),
            decode2IBack = head(decoded$disp.info$decode2IBack, iterations - 1),
            decode2 = tail(decoded$disp.info$decode2Back, 1)[[1]],
            decode2I = tail(decoded$disp.info$decode2IBack, 1)[[1]],
            tempResultSoft = decoded$disp.info$tempResultSoft,
            tempResultHard = decoded$disp.info$tempResultHard,
            result = message.decoded,
            iterations = iterations,
            permutation = permutation.vector,
            encoder = coder.info),
          encoding = "UTF-8")
        rstudioapi::viewer(system.file("pdf", "TurboDecode.pdf", package = "channelcoding"))
      }
    }

    return(message.decoded)
  }

#' Function to generate a permutation vector.
#'
#' This function is a helper function which helps the user to create different permutation
#' vectors. With the argument \code{type} the user can choose one of 6 different interleaver types.
#' Some interleaver need additional informations in the \code{args} argument.
#'
#'
#' \itemize{
#'   \item RANDOM: Creates a random permutation vector. No arguments in \code{args}.
#'   \item PRIMITIVE: Shift the initial vector (\code{c(1,2,3,...)}) so that the \code{args$root}
#'                    value is the index of the 1.
#'   \item CYCLIC: Creates a initial matrix with the arguments \code{args$cols}, \code{args$rows}
#'                 and shift each row of the matrix by the index of the row multiplied by the \code{args$distance}.
#'                 The resulting permutation vector is read out from top to bottom from each column of the matrix.
#'   \item BLOCK: This is a type of interleaver in which the bits are read in from left to right in each row from the matrix
#'                and read out from top to bottom each column. With the arguments \code{args$cols}, \code{args$rows}
#'                the user can change the matrix size.
#'   \item HELICAL: This is a type of interleaver in which the bits are read in from left to right in each row from the matrix
#'                and read out from left top to bottom right. When the last row of the matrix is arrived
#'                then the next bit is read out from the first row but from the next column. With the arguments \code{args$cols}, \code{args$rows}
#'                the user can change the matrix size.
#'   \item DIAGNOAL: The difference to the HELICAL interleaver is that when the last row is arrived
#'                   the next bit is in the first row and the first unused column. With the arguments \code{args$cols}, \code{args$rows}
#'                the user can change the matrix size.
#' }
#'
#' @param message.length Length of message which will be encoded.
#' @param coder.info Coder which will be used for encoding and decoding.
#' @param type Type of the interleaver, possibilities: RANDOM, PRIMITIVE, CYCLIC, BLOCK
#'             HELICAL, DIAGONAL.
#' @param args Arguments for some interleaver. (must be a list)
#' @param visualize Flag to visualize the resulting permutation matrix/vector.
#'
#' @return Created permutation vector.
#'
#' @examples
#' input <- c(1,0,1,1,0,1)
#' coder <- ConvGenerateRscEncoder(2, 2, c(5, 7))
#'
#' #RANDOM
#' permutation <- TurboGetPermutation(length(input), coder, "RANDOM")
#'
#' #PRIMITIV
#' permutation <- TurboGetPermutation(length(input), coder, "PRIMITIVE", list(root=2))
#'
#' #CYCLIC
#' permutation <- TurboGetPermutation(length(input), coder, "CYCLIC", list(cols=2, rows=4, distance=2))
#'
#' #BLOCK
#' permutation <- TurboGetPermutation(length(input), coder, "BLOCK", list(cols=2, rows=4))
#'
#' #HELICAL
#' permutation <- TurboGetPermutation(length(input), coder, "HELICAL", list(cols=2, rows=4))
#'
#' #DIAGONAL
#' permutation <- TurboGetPermutation(length(input), coder, "DIAGONAL", list(cols=2, rows=4))
#'
#' @export
TurboGetPermutation <- function(message.length,
                                coder.info,
                                type = "RANDOM",
                                args = NULL,
                                visualize = FALSE) {

  stopifnot(message.length > 0)

  CheckCoder(coder.info)

  if (!is.character(type)) {
    stop("Type must be an string!")
  }

  switch(type,
    RANDOM = {
      interleaver <- sample(c(0:(message.length + coder.info$M - 1)))

      if (visualize) {
        print(paste("Interleaver-Vector: ", type))
        print(interleaver)
      }

      return(interleaver)
    },
    PRIMITIVE = {
      if (is.null(args$root)) {
        stop("Argument args$root was not set!")
      }
      N <- message.length + coder.info$M - 1
      init <- c(0:N)
      interleaver <- (init - args$root) %% (N + 1)

      if (visualize) {
        print(paste("Interleaver-Vector: ", type))
        print(interleaver)
      }

      return(interleaver)
    },
    CYCLIC = {
      if (is.null(args$rows) |
          is.null(args$cols) | is.null(args$distance)) {
        stop("Argumente (args) was not set correctly!")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (message.length + coder.info$M)) {
        stop("Length of input does not match the amount of rows and columns!")
      }
      N <- rows * cols
      init <- matrix(c(0:(N - 1)), nrow = rows, byrow = FALSE)

      i <- 0
      interleaver <-
        t(apply(init, 1,
          function(x) {
            temp <- Shift(x, args$distance * (i))
            i <<- i + 1
            return(temp)
            }
          ))

      if (visualize) {
        print("Initial-Matrix")
        print(init)
        print(paste("Interleaver-Matrix: ", type))
        print(interleaver)
      }

      return(as.vector(interleaver))
    },
    BLOCK = {
      if (is.null(args$rows) | is.null(args$cols)) {
        stop("Argumente (args) was not set correctly!")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (message.length + coder.info$M)) {
        stop("Length of input does not match the amount of rows and columns!")
      }
      N <- rows * cols
      init <- matrix(c(0:(N - 1)), nrow = rows, byrow = TRUE)

      if (visualize) {
        print("Initial-Matrix")
        print(init)
        print(paste("Interleaver-Vector: ", type))
        print(as.vector((init)))
      }

      return(as.vector((init)))
    },
    HELICAL = {
      if (is.null(args$rows) | is.null(args$cols)) {
        stop("Argumente (args) was not set correctly!")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (message.length + coder.info$M)) {
        stop("Length of input does not match the amount of rows and columns!")
      }
      N <- rows * cols
      init <- vector(mode = "numeric", length = N)

      i <- 0
      interleaver <-
        sapply(init, function(x) {
          x <- (((i %% cols) + (i * cols)) %% N)
          i <<- i + 1
          return(x)
        })

      if (visualize) {
        print("Initial-Matrix")
        print(matrix(c(0:(N - 1)), nrow = rows, byrow = TRUE))
        print(paste("Interleaver-Vector: ", type))
        print(interleaver)
      }

      return(interleaver)
    },
    DIAGONAL = {
      if (is.null(args$rows) | is.null(args$cols)) {
        stop("Argumente (args) was not set correctly!")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (message.length + coder.info$M)) {
        stop("Length of input does not match the amount of rows and columns!")
      }
      N <- rows * cols
      init <- vector(mode = "numeric", length = N)

      i <- 0
      interleaver <-
        sapply(init, function(x) {
          x <- (i * cols) %% N + (i %/% rows + i %% rows) %% cols
          i <<- i + 1
          return(x)
        })

      if (visualize) {
        print("Initial-Matrix")
        print(matrix(c(0:(N - 1)), nrow = rows, byrow = TRUE))
        print(paste("Interleaver-Vector: ", type))
        print(interleaver)
      }

      return(interleaver)
    }
  )
  stop("Type of interleaver (type) was not set correctly!")
}

#' Create a punctuation matrix from a vector.
#'
#' This is a helper function to create a correct punctuation matrix. The punctuation.vector
#' should only contain 0 and 1 and has to be a length which is a multiple of 3. The resulting
#' punctuation matrix should not contain a 0 column because this could cause an error during
#' decoding.
#'
#' @param punctuation.vector Vector which will be mapped to the matrix.
#' @param visualize Flag to visualize the resulting punctuation matrix.
#'
#' @return Punctuation matrix for encoder and decoder.
#'
#' @examples
#' TurboGetPunctuationMatrix(c(1, 1, 0, 1, 0, 1), TRUE)
#'
#' @export
TurboGetPunctuationMatrix <- function(punctuation.vector, visualize = FALSE) {
  if (any((punctuation.vector != 1)[punctuation.vector != 0])) {
    stop("Invalid punctuation vector! Only values 0/1 are allowed!")
  }
  if (length(punctuation.vector) %% 3 != 0) {
    stop("Wrong length of punctuation vector! Must be a multiple of 3!")
  }

  mat <- matrix(punctuation.vector, nrow = 3)

  if (visualize) {
    print("Punctuation-Matrix:")
    print(mat)
  }

  if (any(colSums(mat) == 0)) {
    stop("Punctuation matrix should not have a 0 column!")
  }

  return(mat)
}

#' Function to make a automatic simulation of an coder with different signal/noise ratio.
#'
#' This easy functions makes it possible to compare different coders. The function encode
#' a random message, apply noise to the code and then decode the code. This will be computed
#' many times and after all iterations the bit error rate will be calculated. This procedure
#' is applied to different signal/noise ratios. The result will be printed in a graph, when
#' visualization flag is set to TRUE.
#'
#' @param coder Coder which will be created with \code{\link{ConvGenerateEncoder}} or \code{\link{ConvGenerateRscEncoder}}.
#' @param permutation.type Type of permutation vector.
#' @param permutation.args Arguments to the \code{\link{TurboGetPermutation}} function.
#' @param decode.iterations Amount of decoding iterations inside the turbo decoder.
#' @param msg.length Length of the randomly created message.
#' @param min.db Start value of the signal/noise ratio.
#' @param max.db End value of the signal/noise ration.
#' @param db.interval Interval which will be added to the actual signal/noise ratio after all
#'                    iterations are applied.
#' @param iterations.per.db Amount of iterations each signal/noise ration step.
#' @param punctuation.matrix Punctuation matrix to puncture the output, will be created with \code{\link{TurboGetPunctuationMatrix}}.
#' @param visualize Flag to decide whether to create a visualization pdf or not.
#'
#' @return DataFrame which contains the bit error rate for each signal/noise ratio step.
#'
#' @examples
#' #all default parameters
#' TurboSimulation()
#'
#' #without punctuation
#' coder <- ConvGenerateRscEncoder(2, 2, c(5, 7))
#' TurboSimulation(coder, "RANDOM", NULL, 5, 10, 0.01, 1, 0.05, 50, NULL, FALSE)
#'
#' @export
TurboSimulation <- function(coder = NULL,
                            permutation.type = "PRIMITIVE",
                            permutation.args = list(root=0),
                            decode.iterations = 5,
                            msg.length = 100,
                            min.db = 0.1,
                            max.db = 2.0,
                            db.interval = 0.1,
                            iterations.per.db = 100,
                            punctuation.matrix = NULL,
                            visualize = FALSE)
{
  stopifnot(decode.iterations > 0, msg.length > 0, iterations.per.db > 0,
            min.db > 0, max.db > 0, max.db >= min.db, db.interval > 0)

  if (is.null(coder)) {
    warning("Standard-Coder was used! RSC, N=2, M=2, Generators: (5,7)")
    coder <- ConvGenerateRscEncoder(2,2,c(5,7))
  }

  v.db <- seq(from = min.db, to = max.db, by = db.interval)
  v.ber <- numeric(0)

  perm <- TurboGetPermutation(msg.length, coder, permutation.type, permutation.args)

  total.errors <- 0

  total.iterations <- length(v.db) * iterations.per.db

  print("Turbo Simulation")
  progress.bar <- txtProgressBar(min = 0, max = total.iterations, style = 3)

  progress.counter <- 0

  for (db in v.db) {
    for (i in 1 : iterations.per.db) {
      # create message
      message <- sample(c(0,1), msg.length, replace = TRUE)

      # encode
      coded <- TurboEncode(message, perm, coder, punctuation.matrix = punctuation.matrix)

      # if punctured then take punctured code
      if (!is.null(punctuation.matrix)) {
        coded <- coded$punctured
      }

      # add noise
      noisy <- ApplyNoise(coded, db)

      # decode
      decoded <- TurboDecode(noisy, perm, decode.iterations, coder,
                             punctuation.matrix = punctuation.matrix)

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
    rmarkdown::render(
      system.file("rmd", "Simulation.Rmd", package = "channelcoding"),
      output_dir = system.file("pdf", package = "channelcoding"),
      output_file = "SimulationTurbo.pdf",
      params = list(
        turbo = TRUE,
        message.length = msg.length,
        iterations.per.db = iterations.per.db,
        min.db = min.db,
        max.db = max.db,
        db.interval = db.interval,
        permutation = perm,
        decode.iterations = decode.iterations,
        punctuation = punctuation.matrix,
        encoder = coder,
        dataframe = df),
      encoding = "UTF-8")
    rstudioapi::viewer(system.file("pdf", "SimulationTurbo.pdf", package = "channelcoding"))
  }

  return(df)
}

#' Open the beforehand created visualization PDF files again.
#'
#' With this function it is easy to reopen the PDF files which will be created with
#' \code{\link{TurboEncode}} and \code{\link{TurboDecode}}. The files are stored in the
#' program files of R. (example path: "C:/Program Files/R/R-3.2.4/library/channelcoding/pdf")
#'
#' @param encode Flag to open the encode pdfs.
#' @param punctured Flag to open the decode pdfs.
#' @param simulation Flag to open the simulation pdf.
#'
#' @examples
#' # open encode without punctuation PDF
#' TurboOpenPDF()
#'
#' # open encode with punctuation PDF
#' TurboOpenPDF(punctured = FALSE)
#'
#' # open decode with punctuation PDF
#' TurboOpenPDF(encode = FALSE, punctured = TRUE)
#'
#' # open decode without punctuation PDF
#' TurboOpenPDF(encode = FALSE)
#'
#' # open simulation PDF
#' TurboOpenPDF(simulation = TRUE)
#'
#' @export
TurboOpenPDF <- function(encode = TRUE, punctured = FALSE, simulation = FALSE) {
  if (simulation) {
    path <- system.file("pdf", "SimulationTurbo.pdf", package = "channelcoding")
  } else {
    if (encode) {
      if (punctured) {
        path <- system.file("pdf", "TurboEncodePunctured.pdf", package = "channelcoding")
      } else {
        path <- system.file("pdf", "TurboEncode.pdf", package = "channelcoding")
      }
    } else {
      if (punctured) {
        path <- system.file("pdf", "TurboDecodePunctured.pdf", package = "channelcoding")
      } else {
        path <- system.file("pdf", "TurboDecode.pdf", package = "channelcoding")
      }
    }
  }
  if (path != "") {
    rstudioapi::viewer(path)
  } else {
    warning("File does not exists!")
  }
}
