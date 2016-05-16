#' Message distortion.
#'
#' This function will randomly distort a clean message.
#' @author Bene Wimmer
#' @param msg Message to alter.
#' @param SNR.db Signal-Noise-Ratio in dB simulating the noisy channel.
#' @param binary False = Soft Output, True = Hard Output
#' @return Distorted message containing noise.
#' @export
ApplyNoise <- function(msg, SNR.db = 3, binary = FALSE)
{
  msg.len <- length(msg);
  SNR.linear <- 10^(SNR.db/10);
  power <- sum(msg^2)/(msg.len); #power of vector msg

  # noise is a vector of lenth msg_len and contains
  # normal distributed values with mean 0 and standard-deviation 1
  noise <- sqrt(power / SNR.linear) * rnorm(msg.len,0,1)

  msg.out <- msg + noise;

  if(binary){
    msg.out <- ifelse(msg.out <= 0.5, 0, 1)
  }

  return(msg.out);
}

