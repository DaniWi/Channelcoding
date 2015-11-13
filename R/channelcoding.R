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
  return(msg)
}

