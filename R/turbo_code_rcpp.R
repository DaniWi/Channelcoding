# turbo_code_rcpp.R
#
# interface for convolutional codes
# provided functions:
#  - generate nsc coder
#  - generate rsc coder
#  - encode
#  - decode (soft in & out)
#  - decode (hard decision)

#library(Rcpp)
#sourceCpp('src/turbo_code.cpp');


#' Turbo Encode
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
#' @export
#' @useDynLib channelcoding
TurboEncode <-
  function(message, permutation, encoder.info, parity.index = encoder.info$N) {
    if (encoder.info$N < 2) {
      stop("Error: Encoder muss 2 Ausgängen besitzen!")
    }
    if (encoder.info$generators[1] != 1) {
      stop(
        "Error: Encoder muss ein systematischer Encoder sein! (Ausgang 1 muss Polynom 1 besitzen)"
      )
    }
    if (length(permutation) != (length(message) + encoder.info$M)) {
      stop("Error: Permutation hat die falsche Länge!")
    }

    parity.1 <- conv_encode(message, encoder.info, TRUE)
    message.perm <- message[permutation + 1]
    parity.2 <- conv_encode(message.perm, encoder.info, FALSE)

    if (parity.index > encoder.info$N) {
      parity.index = encoder.info$N
    }

    temp.index <-
      c(rep(FALSE,0), TRUE, rep(FALSE, encoder.info$N - 1))
    message.encoded <- parity.1[temp.index]

    temp.index <-
      c(rep(FALSE,parity.index - 1), TRUE, rep(FALSE, encoder.info$N - parity.index))
    parity.1 <- parity.1[temp.index]
    parity.2 <- parity.2[temp.index]

    return(c(message.encoded, parity.1, parity.2))
  }


#' @export
#' @useDynLib channelcoding
TurboDecode <-
  function(message, permutation, iterations, encoder.info, parity.index) {
    if ((length(message) %% 3) != 0) {
      stop("Error: Nachricht hat die falsche Länge")
    }
    if (length(permutation) != (length(message) / 3)) {
      stop("Error: Permutation hat die falsche Länge!")
    }
    if (encoder.info$generators[1] != 1) {
      stop(
        "Error: Encoder muss ein systematischer Encoder sein! (Ausgang 1 muss Polynom 1 besitzen)"
      )
    }
    if (encoder.info$N < 2) {
      stop("Error: Encoder muss 2 Ausgängen besitzen!")
    }

    message.length <- length(message) / 3
    message.orig <- message[1:message.length]
    parity.1 <- message[(message.length + 1):(2 * message.length)]
    parity.2 <-
      message[(2 * message.length + 1):(3 * message.length)]

    message.decoded.term <-
      c_turbo_decode(
        message.orig, parity.1, parity.2, permutation,
        iterations, encoder.info$N, encoder.info$M, encoder.info$prevState,
        encoder.info$output, parity.index
      )

    output.soft <-
      message.decoded.term$softOutput[1:(message.length - encoder.info$M)]
    output.hard <-
      message.decoded.term$hardOutput[1:(message.length - encoder.info$M)]
    message.decoded <-
      list(output.soft = output.soft, output.hard = output.hard)

    return(message.decoded)
  }








#' @export
#' @useDynLib channelcoding
TurboGetPermutation <- function(length, encoder_info, type, args) {
  switch(type,
         RANDOM = {
           if (is.null(encoder_info$M)) {
             stop("Error: Encoder nicht richtig gesetzt!")
           }
           return(sample(c(0:(
             length + encoder_info$M - 1
           ))))
         },
         PRIMITIVE = {
           if (is.null(args$root)) {
             stop("Error: root not set")
           }
           N <- length + encoder_info$M - 1
           init <- c(0:N)
           interleaver <- (init - args$root) %% (N + 1)
           return(interleaver)
         },
         CYCLIC = {
           if (is.null(args$rows) |
               is.null(args$cols) | is.null(args$distance)) {
             stop("Error: args not set")
           }
           rows <- args$rows
           cols <- args$cols
           if (rows * cols != (length + encoder_info$M)) {
             stop("Error: length of input not correct with rows and cols")
           }
           N <- rows * cols
           init <- matrix(c(0:(N - 1)), nrow = rows, byrow = FALSE)
           print("Original")
           print(init)
           i <- 0
           print("Interleaver Matrix")
           interleaver <-
             t(apply(init,1,function(x) {
               temp <-
                 shift(x,(-1) * args$distance * (i)); i <<- i + 1; return(temp)
             }))
           print(interleaver)
           return(as.vector(interleaver))
         })
}
