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
Encode = function(msg, type, params, visualize)
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
Decode = function(msg, type, params, visualize)
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


#' Message distortion.
#'
#' This function will randomly distort a clean message.
#' @author Bene Wimmer
#' @param msg Message to alter.
#' @param SNR.db Signal-Noise-Ratio in dB simulating the noisy channel.
#' @return Distorted message containing noise.
#' @export
ApplyNoise <- function(msg, SNR.db = 3)
{
  msg.len <- length(msg);
  SNR.linear <- 10^(SNR.db/10);
  power <- sum(msg^2)/(msg.len); #power of vector msg

  # noise is a vector of lenth msg_len and contains
  # normal distributed values with mean 0 and standard-deviation 1
  noise <- sqrt(power / SNR.linear) * rnorm(msg.len,0,1)

  msg.out <- msg + noise;

  return(msg.out);
}

#' Testing BCH codes
#'
#'
#' @author Bene Wimmer
#' @param data MSG
#' @param length Length of MSG
#' @param t Errors to correct
#' @return String
#' @export
testBCH = function(data, length, t, visualize=FALSE)
{
  if(!visualize){
    data = rawToBits(charToRaw(data))
  }

  enc = .BCHEncode(as.integer(data),c(length,t),FALSE)
  enc_error = ifelse(ApplyNoise(enc) <= 0.5, 0, 1)
  dec = .BCHDecode(enc_error,c(length,t),FALSE)

  if(visualize){
    rmarkdown::render(system.file("rmd", "BlockEncode.Rmd", package = "channelcoding"),
                      output_dir = system.file("pdf", package = "channelcoding"),
                      encoding = "UTF-8",
                      params = list(message = data,
                                    enc = enc,
                                    enc_error = enc_error,
                                    dec = dec))

    rstudioapi::viewer(system.file("pdf", "BlockEncode.pdf", package = "channelcoding"))
    return(dec)
  }else{
    return(rawToChar(packBits(dec[1:(length(dec)- (length(dec)%%8))])))
  }



}
