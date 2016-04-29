#' generic encode function
#'
#' Is used as a generic encode funtion that will forward the call depending on the type
#' @author Bene Wimmer
#' @param msg message to encode
#' @param type the type of code you want to use
#' @param params A vector with arguments
#' @param visualize A flag for enabling visualization
#' @return encoded message
#' @export
encode = function(msg, type, params, visualize)
{
  switch(type,
         HAMMING={
           return(.HammingEncode(msg,params,visualize))
         },
         MATRIX={
           return(.MatrixEncode(msg,params,visualize))
         },
         {
           print("No or Wrong Type selected")
           return()
         })

}


#' general decode function
#'
#' Is used as a generic decode funtion that will forward the call depending on the type
#' @author Bene Wimmer
#' @param msg message to decode
#' @param type the type of code you want to use
#' @param params A vector with arguments
#' @param visualize A flag for enabling visualization
#' @return decoded message
#' @export
decode = function(msg, type, params, visualize)
{
  switch(type,
         HAMMING={
           return(.HammingDecode(msg,params,visualize))
         },
         MATRIX={
           return(.MatrixDecode(msg,params,visualize))
         },
         {
           print("No or Wrong Type selected")
           return()
         })
}


#' Failure function
#'
#' This function will randomly disturb a clean message
#' @author Bene Wimmer
#' @param msg message to alter
#' @param SNR.db signal noise ratio of the simulated channel in dB
#' @param visualize A flag for enabling visualization
#' @return message with errors
#' @export
applyNoise <- function(msg, SNR.db = 5, visualize = FALSE)
{
   msg.len <- length(msg);
   SNR.linear <- 10^(SNR.db/10);
   power <- sum(msg^2)/(msg.len); #power of vector msg

   # noise is a vector of lenth msg_len and contains
   # normal distributed values with mean 0 and standard-deviation 1
   noise <- sqrt(power / SNR.linear) * rnorm(msg.len,0,1)

   msg.out <- msg + noise;

   if (visualize) {
      v.flipped <- (msg + msg.out) %% 2;
      # rle: run-length-encoding of the noise(v_flipped) vector
      # results in a table with values and lengths
      rlenc <- rle(v.flipped);
      # get number of n-bit errors (where value==1)
      my.sample <- rlenc$lengths[rlenc$values == 1];
      # plot as table of n-bit error occurences
      # from 1 to highest occuring error (most neighbouring error bits)
      error.sum <- sum(v.flipped);
      max.errorbits <- max(3,rlenc$lengths[rlenc$values == 1]);
      barplot(table(factor(my.sample,levels=1:max.errorbits)),main=paste("Total error bits:",error.sum));
   }

   return(msg.out);
}

TurboSimulation <- function(coder,
                            permutation.type = "RANDOM",
                            decode.iterations = 10,
                            msg.length = 100,
                            iterations.per.db = 100,
                            min.db = 0.1,
                            max.db = 5.0,
                            db.interval = 0.1,
                            punctuation.matrix = NULL)
{
  stopifnot(decode.iterations > 0, msg.length > 0, iterations.per.db > 0,
            min.db > 0, max.db > 0, max.db >= min.db, db.interval > 0)

  v.db <- seq(from = min.db, to = max.db, by = db.interval)

  total.errors <- 0

  for (db in v.db) {
    for (i in 1 : iterations.per.db) {
      # message erzeugen
      message <- sample(c(0,1), msg.length, replace = TRUE)

      # encode
      perm <- TurboGetPermutation(msg.length, encoder, permutation.type)
      coded <- TurboEncode(message, perm, coder, punctuation.matrix = punctuation.matrix)

      # noise hinzufügen
      noisy <- applyNoise(coded, db)

      # anzahl flipped bits (channel errors)
      coded.hard <- ifelse(coded >= 0, 0, 1)
      noisy.hard <- ifelse(noisy >= 0, 0, 1)
      channel.errors <- sum(abs(coded.hard - noisy.hard))

      # decode
      decoded <- TurboDecode(noisy, perm, decode.iterations, coder,
                             punctuation.matrix = punctuation.matrix)

      # vgl decoded & message
      decode.erros <- sum(abs(decoded - message))
    }
  }


  df <- data.frame(db = v.db, ber = v.ber)

}
