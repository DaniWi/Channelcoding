# turbo_code_rcpp.R
#
# interface for turbo codes
# provided functions:
#  - TurboEncode
#  - TurboDecode
#  - TurboGetPermutation


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
#' @param permutation Permutationsvektor der mittels \code{\link{TurboGetPermutation}} erzeugt werden sollte
#' @param coder.info Kodierer
#' @param parity.index Index des zu verwendenten Ausgangs des Kodieres
#'
#' @return Kodierte Nachricht, die aus Ausgangsnachricht und den beiden kodierten Nachrichten besteht
#'
#' @examples
#' input <- c(1,0,1,1,0)
#' coder <- GenerateRscEncoder(2,2,c(5,7))
#' perm <- TurboGetPermutation(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#'
#' @export
TurboEncode <-
  function(message, permutation, coder.info, parity.index = coder.info$N, dotting.matrix = NULL) {
    #Checks all parameters
    if (any((message != 1)[message != 0])) {
      stop("Error: Nachricht darf nur 0er und 1er besitzen!")
    }
    if (coder.info$N < 2) {
      stop("Error: coder muss 2 Ausgängen besitzen!")
    }
    if (coder.info$rsc != TRUE && coder.info$generators[1] != (2^coder.info$M)) {
      stop(
        "Error: coder muss ein systematischer coder sein! (Ausgang 1 muss Polynom 2^M besitzen)"
      )
    }
    if (length(permutation) != (length(message) + coder.info$M)) {
      stop("Error: Permutation hat die falsche Länge!")
    }
    if(!is.null(dotting.matrix) && ((dim(dotting.matrix))[1] != 3)) {
      stop("Error: Punktierungsmatrix hat falsche Anzahl an Zeilen, bei Turbo-Codes müssen es 3 sein!")
    }

    #when parity.index is too high
    if (parity.index > coder.info$N) {
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
    if(!is.null(dotting.matrix)) {
      #dotting the output
      result.dot <- PunctureCode(result.orig, dotting.matrix)
    }

    #rmarkdown::render(system.file("rmd", "TurboEncode.Rmd", package = "channelcoding"), params = list(
    #  orig = message,
    #  interl = message.perm,
    #  parity1 = parity.1,
    #  parity2 = parity.2,
    #  result = c(message.encoded, parity.1, parity.2)))
    #rstudioapi::viewer(system.file("rmd", "TurboEncode.pdf", package = "channelcoding"))

    if(!is.null(dotting.matrix)) {
      return(list(original=result.orig, punctured=result.dot))
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
#' @param permutation Permutationsvektor der mittels \code{\link{TurboGetPermutation}} erzeugt werden sollte
#' @param coder.info Dekodierer
#' @param parity.index Index des zu verwendenten Ausgangs des Kodieres
#'
#' @return Deodierte Nachricht
#'
#' @examples
#' input <- c(1,0,1,1,0)
#' coder <- GenerateRsccoder(2,2,c(5,7))
#' perm <- TurboGetPermutation(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#' result <- TurboDecode(input, perm, 5, coder)
#'
#' @export
TurboDecode <-
  function(message, permutation, iterations, coder.info, parity.index = coder.info$N, dotting.matrix = NULL) {
    #Checks all parameters
    if (coder.info$N < 2) {
      stop("Error: coder muss 2 Ausgängen besitzen!")
    }
    if (coder.info$rsc != TRUE && coder.info$generators[1] != (2^coder.info$M)) {
      stop(
        "Error: coder muss ein systematischer coder sein! (Ausgang 1 muss Polynom 2^M besitzen)"
      )
    }
    if(!is.null(dotting.matrix) && ((dim(dotting.matrix))[1] != 3)) {
      stop("Error: Punktierungsmatrix hat falsche Anzahl an Zeilen, bei Turbo-Codes müssen es 3 sein!")
    }

    if (parity.index > coder.info$N) {
      parity.index = coder.info$N
    }


    if(!is.null(dotting.matrix)) {
      #insert missing bits from puncturing
      message.afterPunct <- InsertPuncturingBits(message, dotting.matrix)
    } else {
      message.afterPunct <- message
    }

    #Check length of input(with inserted bits) and permutation
    if ((length(message.afterPunct) %% 3) != 0) {
      if(!is.null(dotting.matrix)) {
        stop("Error: Punktierung hat fehlgeschlagen, andere Punktierungsmatrix wählen")
      }
      stop("Error: Nachricht hat die falsche Länge")
    }
    if (length(permutation) != (length(message.afterPunct) / 3)) {
      stop("Error: Permutation hat die falsche Länge!")
    }

    #parse original message from input message
    message.length <- length(message.afterPunct) / 3
    message.orig <- Deinterleave(message.afterPunct, 1)
    parity.1 <- Deinterleave(message.afterPunct, 2)
    parity.2 <- Deinterleave(message.afterPunct, 3)

    #decode message
    message.decoded.term <-
      c_turbo_decode(
        message.orig, parity.1, parity.2, permutation,
        iterations, coder.info$N, coder.info$M, coder.info$prev.state,
        coder.info$output, parity.index
      )

    #delete termination bits from result
    output.soft <-
      message.decoded.term$soft.output[1:(message.length - coder.info$M)]
    output.hard <-
      message.decoded.term$hard.output[1:(message.length - coder.info$M)]
    message.decoded <-
      list(output.soft = output.soft, output.hard = output.hard)

    return(message.decoded)
  }

#' @export
TurboGetPermutation <- function(length, coder.info, type, args) {
  if (is.null(coder.info$M)) {
    stop("Error: Kodierer nicht richtig gesetzt!")
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
        stop("Error: root(args) wurde nicht gesetzt")
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
        stop("Error: Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Error: Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
      }
      N <- rows * cols
      init <- matrix(c(0:(N - 1)), nrow = rows, byrow = FALSE)

      i <- 0
      interleaver <-
        t(apply(init,1,function(x) {
          temp <-
            binhf::shift(x,args$distance * (i)); i <<-
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
        stop("Error: Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Error: Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
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
        stop("Error: Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Error: Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
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
        stop("Error: Argumente wurden nicht richtig gesetzt")
      }
      rows <- args$rows
      cols <- args$cols
      if (rows * cols != (length + coder.info$M)) {
        stop("Error: Länge von Input stimmt nicht mit Reihen- und Spaltenanzahl zusammen!")
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
  stop("Error: Type von Interleaver wurde nicht richtig gewählt!")
}

#' @export
TurboGetPuncturingMatrix <- function(puncturing.vector) {

  if (any((puncturing.vector != 1)[puncturing.vector != 0])) {
    # puncturing.vector has elements with value different from 0 or 1 which is not allowed
    stop("Invalid puncturing vector! Only values 0/1 are allowed!")
  }

  return(matrix(puncturing.vector, nrow = 3))
}
