#use of .C
dyn.load("src/helloC.dll")
helloC <- function(greeting) {
  if (!is.character(greeting)) {
    stop("Argument 'greeting' must be of type 'character'.")
  }
  result <- .C("helloC",
               greeting=greeting,
               count=as.integer(1))
  return(result$count)
}

#use of .Call
dyn.load("src/helloD.dll")
helloD <- function(greeting) {
  result <- .Call("helloD", greeting)
  return(result)
}

#use of .Call with cpp file
Rcpp::sourceCpp('src/helloE.cpp')
helloE <- function(greeting) {
  result <- helloEcpp(greeting)
}
