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
         BCH={
           return(.BCHEncode(msg,params,visualize))
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
         BCH={
           return(.BCHDecode(msg,params,visualize))
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
#' @param params A vector with arguments
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

#' Testing BCH codes
#'
#'
#' @author Bene Wimmer
#' @param String
#' @param length
#' @param t
#' @return String
#' @export
testBCH = function(data, length, t)
{
  enc = .BCHEncode(as.integer(rawToBits(charToRaw(data))),c(length,t),FALSE)

  dec = .BCHDecode(applyNoise(enc),c(length,t),FALSE)


  return(rawToChar(packBits(dec[1:(length(dec)- (length(dec)%%8))])))

}
