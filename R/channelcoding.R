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
#' @param params A vector with arguments
#' @param visualize A flag for enabling visualization
#' @return message with errors
#' @export
applyNoise = function(msg, params, visualize = FALSE)
{

   SNR_db = params; # Signal-Noise-Ratio in dB
   msg_len = length(msg);
   SNR_linear = 10^(SNR_db/10);
   power = sum(msg^2)/(msg_len); #power of vector msg

   # noise is a vector of lenth msg_len and contains
   # normal distributed values with mean 0 and standard-deviation 1
   noise = sqrt(power / SNR_linear) * rnorm(msg_len,0,1)

   msg_out = msg + noise;

   # map every vector element <= 0.5 to 0 and > 0.5 to 1
   #msg_out = ifelse(msg_out <= 0.5, 0, 1)

   if (visualize) {
      v_flipped <- (msg + msg_out) %% 2;
      # rle: run-length-encoding of the noise(v_flipped) vector
      # results in a table with values and lengths
      rlenc <- rle(v_flipped);
      # get number of n-bit errors (where value==1)
      mysample <- rlenc$lengths[rlenc$values == 1];
      # plot as table of n-bit error occurences
      # from 1 to highest occuring error (most neighbouring error bits)
      error_sum <- sum(v_flipped);
      max_errorbits <- max(3,rlenc$lengths[rlenc$values == 1]);
      barplot(table(factor(mysample,levels=1:max_errorbits)),main=paste("Total error bits:",error_sum));
   }

   return(msg_out);
}
