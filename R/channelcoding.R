#' general encdoe funtion
#'
#' @author Witsch Daniel
#' @param msg message to encode
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

decode = function(msg, type, params, visualize)
{
  return(msg)
}

applyNoise = function(msg, params, visualize)
{
   SNR_db = 10;
   msg_len = length(msg);
   
   SNR_linear = 10^(SNR_db/10);
   
   power = sum(msg)/(msg_len); #power of vector msg
   
   noise = sqrt(power / SNR_linear) * rnorm(msg_len,0,1)
   
   msg_out = msg + noise;
   
   
   return(msg_out);
}

