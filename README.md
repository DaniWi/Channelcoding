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
