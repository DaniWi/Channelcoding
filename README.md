# Channelcoding
R-Package for channelcoding (block, convolutional, turbo)

C Code einbinden:
1)C code schreiben
2)R CMD SHLIB file.c ausführen
3)dyn.load("C/helloA.dll") in R ausführen
4)Wrapper Funktion schreiben
	helloB <- function() {
	  result <- .C("helloA",
				   greeting="")
	  return(result$greeting)
	}
5)C Code über Wrapper Funktion ausführen

#' @examples
#' encoder <- generateConvEncoder_nsc(2,2,c(7,5))
#' code <- conv_encode(c(1,0,1), encoder)
#' msg <- conv_decode(code, encoder)

#' @examples
#' encoder <- generateConvEncoder_nsc(2,2,c(7,5))
#' code <- conv_encode(c(1,0,1), encoder)
#' msg <- conv_decode_hard(code, encoder)

#' @examples
#' encoder <- generateConvEncoder_nsc(2,2,c(7,5))
#' code <- conv_encode(c(1,0,1), encoder)

#' @example generateConvEncoder_rsc(2,2,c(1,10,13))
