# turbo_code_rcpp.R
#
# interface for turbo codes
# provided functions:
#  - TurboEncode
#  - TurboDecode
#  - TurboGetpermutation.vector


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
#' perm <- TurboGetpermutation.vector(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#'
#' @export
TurboEncode <-
  function(message, permutation.vector, coder.info, parity.index = coder.info$N, punctuation.matrix = NULL) {
    #Checks all parameters

    stopifnot(length(message) > 0)

    if (any((message != 1)[message != 0])) {
      stop("Nachricht darf nur 0er und 1er enthalten!")
    }
    if (coder.info$N < 2) {
      stop("Kodierer muss mindestens 2 Ausgänge besitzen!")
    }
    if (coder.info$rsc != TRUE && coder.info$generators[1] != (2^coder.info$M)) {
      stop(
        "Kodierer muss ein systematischer coder sein! (Ausgang 1 muss Polynom 2^M besitzen)"
      )
    }
    if (length(permutation.vector) != (length(message) + coder.info$M)) {
      stop("Permutationsvektor hat die falsche Länge!")
    }
    if (any(is.na(match(0:length(permutation.vector)-1)))) {
      stop("Permutationsvektor hat falsche Einträge!")
    }
    if (!is.null(punctuation.matrix) && nrow(punctuation.matrix) != 3) {
      stop("Punktierungsmatrix hat die falsche Anzahl an Zeilen, bei Turbo-Codes müssen es 3 sein!")
    }

    #when parity.index is too high or to low
    if (parity.index > coder.info$N || parity.index < 0) {
      parity.index = coder.info$N
    }

    #first coder with termination
    parity.1 <- ConvEncode(message, coder.info, TRUE)
    #extract original message with termination bits from parity
    temp.index <-
      c(rep(FALSE,0), TRUE, rep(FALSE, coder.info$N - 1))
    message.encoded <- parity.1[temp.index]

    #permutate original message with termination (mapping important!!)
    message.perm <- as.numeric(message.encoded[perm + 1] < 0)

    #second coder without termination
    parity.2 <- ConvEncode(message.perm, coder.info, FALSE)


    #extract parity bits from the output messages
    temp.index <-
      c(rep(FALSE,parity.index - 1), TRUE, rep(FALSE, coder.info$N - parity.index))
    parity.1 <- parity.1[temp.index]
    parity.2 <- parity.2[temp.index]

    result.orig <- Interleave(message.encoded, parity.1, parity.2)
    if(!is.null(punctuation.matrix)) {
      #punctuation the output
      result.punct <- PunctureCode(result.orig, punctuation.matrix)
    }

    #rmarkdown::render(system.file("rmd", "TurboEncode.Rmd", package = "channelcoding"), params = list(
    #  orig = message,
    #  interl = message.perm,
    #  parity1 = parity.1,
    #  parity2 = parity.2,
    #  result = c(message.encoded, parity.1, parity.2)))
    #rstudioapi::viewer(system.file("rmd", "TurboEncode.pdf", package = "channelcoding"))

    if(!is.null(punctuation.matrix)) {
      return(list(original=result.orig, punctured=result.punct))
    } else {
      return(result.orig)
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
#' input <- c(1,0,1,1,0)
#' coder <- GenerateRsccoder(2,2,c(5,7))
#' perm <- TurboGetpermutation.vector(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#' result <- TurboDecode(input, perm, 5, coder)
#'
#' @export
TurboDecode <-
  function(code, permutation.vector, iterations, coder.info, parity.index = coder.info$N, punctuation.matrix = NULL) {
    #Checks all parameters

    stopifnot(length(code) > 0, iterations > 0)

    if (coder.info$N < 2) {
      stop("Kodierer muss 2 Ausgänge besitzen!")
    }
    if (coder.info$rsc != TRUE && coder.info$generators[1] != (2^coder.info$M)) {
      stop(
        "Kodierer muss ein systematischer coder sein! (Ausgang 1 muss Polynom 2^M besitzen)"
      )
    }
    if(!is.null(punctuation.matrix) && nrow(punctuation.matrix) != 3) {
      stop("Punktierungsmatrix hat die falsche Anzahl an Zeilen, bei Turbo-Codes müssen es 3 sein!")
    }

    #when parity.index is too high or to low
    if (parity.index > coder.info$N || parity.index < 0) {
      parity.index = coder.info$N
    }


    if(!is.null(punctuation.matrix)) {
      #insert missing bits from puncturing
      message <- InsertPuncturingBits(code, punctuation.matrix)
    } else {
      message <- code
    }

    #Check length of input(with inserted bits) and permutation.vector
    if ((length(message) %% 3) != 0) {
      if(!is.null(punctuation.matrix)) {
        stop("Fehler während der Punktierung!")
      }
      stop("Code hat die falsche Länge!")
    }
    if (length(permutation.vector) != (length(message) / 3)) {
      stop("permutation.vector hat die falsche Länge!")
    }

    #parse original message from input message
    message.length <- length(message) / 3
    message.orig <- Deinterleave(message, 1)
    parity.1 <- Deinterleave(message, 2)
    parity.2 <- Deinterleave(message, 3)

    #decode message
    message.decoded.term <-
      c_turbo_decode(
        message.orig, parity.1, parity.2, permutation.vector,
        iterations, coder.info$N, coder.info$M, coder.info$prev.state,
        coder.info$output, parity.index
      )

    #delete termination bits from result
    output.soft <-
      head(message.decoded.term$soft.output, message.length - coder.info$M)
    output.hard <-
      head(message.decoded.term$hard.output, message.length - coder.info$M)
    message.decoded <-
      list(output.soft = output.soft, output.hard = output.hard)

    return(message.decoded)
  }

#' @export
TurboGetpermutation.vector <- function(length, coder.info, type, args) {
  if (is.null(coder.info$M)) {
    stop("Kodierer nicht richtig gesetzt!")
  }

  switch(
    type,
    RANDOM = {
      return(sample(c(0:(
        length + coder.info$M - 1
      ))))
    },
    PRIMITIVE = {
      if (is.null(args$root)) {
        stop("root(args) wurde nicht gesetzt")
      }
      N <- length + coder.info$M - 1
      init <- c(0:N)
      interleaver <- (init - args$root) %% (N + 1)

      print("Interleaver Vektor")
      print(interleaver)

      return(interleaver)
    },
    CYCLIC = {
      if (is.null(args$rows) |
          is.null(args$cols) | is.null(args$distance)) {
        stop("Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
      }
      N <- rows * cols
      init <- matrix(c(0:(N - 1)), nrow = rows, byrow = FALSE)

      i <- 0
      interleaver <-
        t(apply(init,1,function(x) {
          temp <-
            Shift(x,args$distance * (i)); i <<-
              i + 1; return(temp)
        }))

      print("Original")
      print(init)
      print("Interleaver Matrix")
      print(interleaver)

      return(as.vector(interleaver))
    },
    BLOCK = {
      if (is.null(args$rows) | is.null(args$cols)) {
        stop("Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
      }
      N <- rows * cols
      init <- matrix(c(0:(N - 1)), nrow = rows, byrow = TRUE)

      print("Original")
      print(init)
      print("Interleaver Vektor")
      print(as.vector((init)))

      return(as.vector((init)))
    },
    HELICAL = {
      if (is.null(args$rows) | is.null(args$cols)) {
        stop("Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
      }
      N <- rows * cols
      init <- vector(mode = "numeric", length = N)

      i <- 0
      interleaver <-
        sapply(init, function(x) {
          x <- (((i %% cols) + (i * cols)) %% 15); i <<- i + 1; return(x)
        })

      print("Original")
      print(matrix(c(0:(N - 1)), nrow = rows, byrow = TRUE))
      print("Interleaver Vektor")
      print(interleaver)

      return(interleaver)
    },
    DIAGONAL = {
      if (is.null(args$rows) | is.null(args$cols)) {
        stop("Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
      }
      N <- rows * cols
      init <- vector(mode = "numeric", length = N)

      i <- 0
      interleaver <-
        sapply(init, function(x) {
          x <- (i * cols) %% N + (i %/% rows + i %% rows) %% cols; i <<- i + 1; return(x)
        })

      print("Original")
      print(matrix(c(0:(N - 1)), nrow = rows, byrow = TRUE))
      print("Interleaver Vektor")
      print(interleaver)

      return(interleaver)
    }
  )
  stop("Type von Interleaver wurde nicht richtig gewählt!")
}

#' @export
TurboGetPuncturingMatrix <- function(puncturing.vector) {

  if (any((puncturing.vector != 1)[puncturing.vector != 0])) {
    # puncturing.vector has elements with value different from 0 or 1 which is not allowed
    stop("Invalid puncturing vector! Only values 0/1 are allowed!")
  }

  return(matrix(puncturing.vector, nrow = 3))
}
