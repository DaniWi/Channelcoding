# turbo_code_rcpp.R
#
# interface for convolutional codes
# provided functions:
#  - generate nsc coder
#  - generate rsc coder
#  - encode
#  - decode (soft in & out)
#  - decode (hard decision)

library(Rcpp)
sourceCpp('src/turbo_code.cpp');

#' @export
#' @useDynLib channelcoding
turbo_encode <- function(message, permutation, encoder_info, parity_index = encoder_info$N) {
  if(encoder_info$N < 2) {
    stop("Error: Encoder muss 2 Ausgängen besitzen!")
  }
  if(encoder_info$generators[1] != 1) {
    stop("Error: Encoder muss ein systematischer Encoder sein! (Ausgang 1 muss Polynom 1 besitzen)")
  }
  if(length(permutation) != (length(message)+encoder_info$M)) {
    stop("Error: Permutation hat die falsche Länge!")
  }

  parity_1 <- conv_encode(message, encoder_info, TRUE)
  message_perm <- message[permutation + 1]
  parity_2 <- conv_encode(message_perm, encoder_info, FALSE)

  if(parity_index > encoder_info$N) {
    parity_index = encoder_info$N
  }

  temp_index <- c(rep(FALSE,0),TRUE,rep(FALSE,encoder_info$N-1))
  message_encoded <- parity_1[temp_index]

  temp_index <- c(rep(FALSE,parity_index-1),TRUE,rep(FALSE,encoder_info$N-parity_index))
  parity_1 <- parity_1[temp_index]
  parity_2 <- parity_2[temp_index]

  return(c(message_encoded,parity_1))
}



#' @export
turbo_get_permutation <- function(length, encoder_info, type, args){

  switch(type,
         RANDOM={
           if(is.null(encoder_info$M)){
             stop("Error: Encoder nicht richtig gesetzt!")
           }
           return(sample(c(0:(length+encoder_info$M-1))))
         },
         PRIMITIVE={
           if(is.null(args$root)){
             stop("Error: root not set")
           }
           N <- length+encoder_info$M-1
           init <- c(0:N)
           interleaver <- (init - args$root) %% (N+1)
           return(interleaver)
         },
         CYCLIC={
           if(is.null(args$rows) | is.null(args$cols) | is.null(args$distance)){
             stop("Error: args not set")
           }
           rows <- args$rows
           cols <- args$cols
           if(rows*cols != (length+encoder_info$M)){
             stop("Error: length of input not correct with rows and cols")
           }
           N <- rows * cols
           init <- matrix(c(0:(N-1)), nrow = rows, byrow = FALSE)
           print("Original")
           print(init)
           i <- 0
           print("Interleaver Matrix")
           interleaver <- t(apply(init,1,function(x) {temp <- shift(x,(-1)*args$distance*(i)); i <<- i + 1; return(temp)}))
           print(interleaver)
           return(as.vector(interleaver))
         })
}
