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
  return(msg)
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
applyNoise = function(msg, params, visualize)
{
   
   SNR_db = 10; # Signal-Noise-Ratio in dB
   msg_len = length(msg);
   SNR_linear = 10^(SNR_db/10);
   power = sum(msg)/(msg_len); #power of vector msg
   
   # noise is a vector of lenth msg_len and contains
   # normal distributed values with mean 0 and standard-deviation 1
   noise = sqrt(power / SNR_linear) * rnorm(msg_len,0,1)
   
   msg_out = msg + noise;
   
   # map every vector element <= 0.5 to 0 and > 0.5 to 1
   msg_out = ifelse(msg_out <= 0.5, 0, 1)
   #msg_out = pmax(pmin(round(msg_out,0),1),0)
   
   return(msg_out);
}

