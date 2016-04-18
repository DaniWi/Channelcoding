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
#' @param encoder.info Kodierer
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
  function(message, permutation, encoder.info, parity.index = encoder.info$N) {
    #Checks all parameters
    if (encoder.info$N < 2) {
      stop("Error: Encoder muss 2 Ausgängen besitzen!")
    }
    if (encoder.info$rsc != TRUE && encoder.info$generators[1] != (2^encoder.info$M)) {
      stop(
        "Error: Encoder muss ein systematischer Encoder sein! (Ausgang 1 muss Polynom 2^M besitzen)"
      )
    }
    if (length(permutation) != (length(message) + encoder.info$M)) {
      stop("Error: Permutation hat die falsche Länge!")
    }

    if (parity.index > encoder.info$N) {
      parity.index = encoder.info$N
    }

    #first encoder with termination
    parity.1 <- ConvEncode(message, encoder.info, TRUE)
    #extract original message with termination bits from parity
    temp.index <-
      c(rep(FALSE,0), TRUE, rep(FALSE, encoder.info$N - 1))
    message.encoded <- parity.1[temp.index]

    #permutate original message with termination
    message.perm <- as.numeric(message.encoded[perm + 1] > 0)

    #second encoder without termination
    parity.2 <- ConvEncode(message.perm, encoder.info, FALSE)


    #extract parity bits from the output messages
    temp.index <-
      c(rep(FALSE,parity.index - 1), TRUE, rep(FALSE, encoder.info$N - parity.index))
    parity.1 <- parity.1[temp.index]
    parity.2 <- parity.2[temp.index]

    rmarkdown::render(system.file("rmd", "TurboEncode.Rmd", package = "channelcoding"), params = list(
      orig = message,
      interl = message.perm,
      parity1 = parity.1,
      parity2 = parity.2,
      result = c(message.encoded, parity.1, parity.2)))
    rstudioapi::viewer(system.file("rmd", "TurboEncode.pdf", package = "channelcoding"))

    return(c(message.encoded, parity.1, parity.2))
  }


#' Dekodieren einer Nachricht mittels dem Turbo-Code-Verfahren
#'
#' Dekodiert eine Nachricht mittels dem Turbo-Code-Verfahren. Dabei wird die Nachricht in 2 systematische Dekodierer
#' gesteckt, wobei bei einem Kodierer die Nachricht permutiert verarbeitet wird. Der mitgegebene Kodierer muss
#' mindestens 2 Ausgänge haben und sollte Systematisch sein. Dadurch wird der erste Ausgang automatisch
#' bei beiden Deodierern durchgeschalten. Nachdem die Nachricht durche beide Dekodierer durch ist, kann dieser Schritt
#' mehrmals wiederholt werden. Mit dem Parameter iterations, kann die Anzahl der Durchläufe verändert werden.
#' Je mehr Durchläufe, desto besser das Ergebnis. Das Prinzip der Turbo-Codes ist, dass die Dekodierung mit Soft-Werten
#' vollzogen wird. Das heißt, dass beim Eingang der tatsächliche Signalpegel berücksichtigt wird und beim Ausgang ein
#' Wahrscheinlichkeitswert berechnet wird, der angibt, wie wahrscheinlich ein Bit am Ausgang ist
#'
#' @author Witsch Daniel
#'
#' @param message Nachricht die dekodiert werden sollte
#' @param permutation Permutationsvektor der mittels \code{\link{TurboGetPermutation}} erzeugt werden sollte
#' @param encoder.info Dekodierer
#' @param parity.index Index des zu verwendenten Ausgangs des Kodieres
#'
#' @return Deodierte Nachricht
#'
#' @examples
#' input <- c(1,0,1,1,0)
#' coder <- GenerateRscEncoder(2,2,c(5,7))
#' perm <- TurboGetPermutation(length(input), coder, "RANDOM")
#' message.encoded <- TurboEncode(input, perm, coder)
#' result <- TurboDecode(input, perm, 5, coder)
#'
#' @export
TurboDecode <-
  function(message, permutation, iterations, encoder.info, parity.index = encoder.info$N) {
    #Checks all parameters
    if (encoder.info$N < 2) {
      stop("Error: Encoder muss 2 Ausgängen besitzen!")
    }
    if (encoder.info$rsc != TRUE && encoder.info$generators[1] != (2^encoder.info$M)) {
      stop(
        "Error: Encoder muss ein systematischer Encoder sein! (Ausgang 1 muss Polynom 2^M besitzen)"
      )
    }
    if ((length(message) %% 3) != 0) {
      stop("Error: Nachricht hat die falsche Länge")
    }
    if (length(permutation) != (length(message) / 3)) {
      stop("Error: Permutation hat die falsche Länge!")
    }

    if (parity.index > encoder.info$N) {
      parity.index = encoder.info$N
    }

    #parse original message from input message
    message.length <- length(message) / 3
    message.orig <- message[1:message.length]

    #parse parity1 from input message
    parity.1 <- message[(message.length + 1):(2 * message.length)]
    #parse parity2 from input message
    parity.2 <-
      message[(2 * message.length + 1):(3 * message.length)]

    #decode message
    message.decoded.term <-
      c_turbo_decode(
        message.orig, parity.1, parity.2, permutation,
        iterations, encoder.info$N, encoder.info$M, encoder.info$prev.state,
        encoder.info$output, parity.index
      )

    #delete termination bits from result
    output.soft <-
      message.decoded.term$soft.output[1:(message.length - encoder.info$M)]
    output.hard <-
      message.decoded.term$hard.output[1:(message.length - encoder.info$M)]
    message.decoded <-
      list(output.soft = output.soft, output.hard = output.hard)

    return(message.decoded)
  }

#' @export
TurboGetPermutation <- function(length, encoder.info, type, args) {
  if (is.null(encoder.info$M)) {
    stop("Error: Encoder nicht richtig gesetzt!")
  }

  switch(
    type,
    RANDOM = {
      return(sample(c(0:(
        length + encoder.info$M - 1
      ))))
    },
    PRIMITIVE = {
      if (is.null(args$root)) {
        stop("Error: root(args) wurde nicht gesetzt")
      }
      N <- length + encoder.info$M - 1
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
      if (rows * cols != (length + encoder.info$M)) {
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
      if (rows * cols != (length + encoder.info$M)) {
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
      if (rows * cols != (length + encoder.info$M)) {
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
      if (rows * cols != (length + encoder.info$M)) {
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
