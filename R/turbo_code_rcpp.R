# turbo_code_rcpp.R
#
# interface for turbo codes
# provided functions:
#  - TurboEncode
#  - TurboDecode
#  - TurboGetpermutation
#  - TurboGetPunctuationMatrix


#' Kodieren einer Nachricht mittels dem Turbo-Code-Verfahren
#'
#' Kodiert eine Nachricht mittels dem Turbo-Code-Verfahren. Dabei wird die Nachricht in 2 systematische Kodierer
#' gesteckt, wobei bei einem Kodierer die Nachricht permutiert verarbeitet wird. Der mitgegebene Kodierer muss
#' mindestens 2 Ausgänge haben und sollte Systematisch sein. Dadurch wird der erste Ausgang automatisch
#' bei beiden Kodierern durchgeschalten, somit muss die Ausgangsnachricht nur einmal übertragen werden,
#' dadurch wird die Koderate verbessert. Um den Minimalabstand der Kodierer zu erhöhen, wird die Nachricht
#' permutiert in den 2ten Kodierer geschickt. Am Ende wird die Ausgangsnachricht mit den beiden kodierten Nachrichten
#' aus den Kodierern verknüpft und retour gegeben. Da nur ein Ausgang von einem Kodierer verwendet wir, muss der Index
#' des Ausgangs angegeben werden.
#'
#' @author Witsch Daniel
#'
#' @param message Nachricht die kodiert werden sollte
#' @param permutation.vector permutation.vectorsvektor der mittels \code{\link{TurboGetpermutation.vector}} erzeugt werden sollte
#' @param coder.info Kodierer
#' @param parity.index Index des zu verwendenten Ausgangs des Kodieres
#'
#' @return Kodierte Nachricht, die aus Ausgangsnachricht und den beiden kodierten Nachrichten besteht
#'
#' @examples
#' input <- c(1,0,1,1,0)
#' coder <- GenerateRscEncoder(2,2,c(5,7))
#' perm <- TurboGetpermutation(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#'
#' punct <- TurboGetPunctuationMatrix(c(1,1,0,1,0,1))
#' message.encoded.punct <- TurboEncode(input, perm, coder, punctuation.matrix = punct)
#'
#' @export
TurboEncode <-
  function(message, permutation.vector = NULL,
           coder.info = NULL,
           parity.index = coder.info$N,
           punctuation.matrix = NULL,
           visualize = FALSE) {

    if (is.null(coder.info)) {
      warning("Standard-Coder was used! RSC, N=2, M=2, Generators: (5,7)")
      coder.info <- GenerateRscEncoder(2,2,c(5,7))
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


#' Dekodieren einer Nachricht mittels dem Turbo-Code-Verfahren
#'
#' Dekodiert eine Nachricht mittels dem Turbo-Code-Verfahren. Dabei wird die Nachricht in 2 systematische Dekodierer
#' gesteckt, wobei bei einem Kodierer die Nachricht permutiert verarbeitet wird. Der mitgegebene Kodierer muss
#' mindestens 2 Ausgänge haben und sollte Systematisch sein. Dadurch wird der erste Ausgang automatisch
#' bei beiden Deodierern durchgeschalten. Nachdem die Nachricht durche beide Dekodierer durch ist, kann dieser Schritt
#' mehrmals wiederholt werden. Mit dem Parameter iterations, kann die Anzahl der Durchläufe verändert werden.
#' Je mehr Durchläufe, desto besser das Ergebnis. Das Prinzip der Turbo-Codes ist, dass die Dekodierung mit Soft-Werten
#' vollzogen wird. Das hei?t, dass beim Eingang der tatsächliche Signalpegel berücksichtigt wird und beim Ausgang ein
#' Wahrscheinlichkeitswert berechnet wird, der angibt, wie wahrscheinlich ein Bit am Ausgang ist
#'
#' @author Witsch Daniel
#'
#' @param message Nachricht die dekodiert werden sollte
#' @param permutation.vector permutation.vectorsvektor der mittels \code{\link{TurboGetpermutation.vector}} erzeugt werden sollte
#' @param coder.info Dekodierer
#' @param parity.index Index des zu verwendenten Ausgangs des Kodieres
#'
#' @return Deodierte Nachricht
#'
#' @examples
#'
#' # Ohne Punktierung
#' input <- c(1,0,1,1,0)
#' coder <- GenerateRsccoder(2,2,c(5,7))
#' perm <- TurboGetpermutation.vector(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#' result <- TurboDecode(message.encoded, perm, 5, coder)
#'
#' # Mit Punktierung
#' punct <- TurboGetPunctuationMatrix(c(1,1,0,1,0,1))
#' message.encoded.punct <- TurboEncode(input, perm, coder, punctuation.matrix = punct)
#' result.punct <- TurboDecode(message.encoded.punct$punctured, perm, 5, coder, punctuation.matrix = punct)
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

    if (is.null(coder.info)) {
      warning("Standard-Coder was used! RSC, N=2, M=2, Generators: (5,7)")
      coder.info <- GenerateRscEncoder(2,2,c(5,7))
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

#' @export
TurboSimulation <- function(coder = NULL,
                            permutation.type = "PRIMITIVE",
                            permutation.args = list(root=0),
                            decode.iterations = 10,
                            msg.length = 100,
                            iterations.per.db = 100,
                            min.db = 0.1,
                            max.db = 2.0,
                            db.interval = 0.1,
                            punctuation.matrix = NULL,
                            visualize = TRUE)
{
  stopifnot(decode.iterations > 0, msg.length > 0, iterations.per.db > 0,
            min.db > 0, max.db > 0, max.db >= min.db, db.interval > 0)

  if (is.null(coder)) {
    warning("Standard-Coder was used! RSC, N=2, M=2, Generators: (5,7)")
    coder <- GenerateRscEncoder(2,2,c(5,7))
  }

  v.db <- seq(from = min.db, to = max.db, by = db.interval)
  v.ber <- numeric(0)

  perm <- TurboGetPermutation(msg.length, coder, permutation.type, permutation.args)

  total.errors <- 0

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
      noisy <- applyNoise(coded, db)

      # decode
      decoded <- TurboDecode(noisy, perm, decode.iterations, coder,
                             punctuation.matrix = punctuation.matrix)

      # vgl decoded & message
      decode.errors <- sum(abs(decoded$output.hard - message))

      total.errors <- total.errors + decode.errors
    }

    v.ber <- c(v.ber, total.errors / (msg.length * iterations.per.db))
    total.errors <- 0
  }

  df <- data.frame(db = v.db, ber = v.ber)

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


  df <- data.frame(db = v.db, ber = v.ber)

  return(df)
}

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
    stop("File does not exists!")
  }
}
