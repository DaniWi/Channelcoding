library(Rcpp);
sourceCpp('src/convolution.cpp');

generateConvEncoder <- function(N, M, generators) {

  matrixList <- c_generateMatrices(N,M,generators);

  convEncoder <- list(N = N,
                      M = M,
                      generators = generators,
                      nextState = matrixList$nextState,
                      prevState = matrixList$prevState,
                      output = matrixList$output);

  return(convEncoder);
}

conv_encode <- function(message, convEncoder) {

  code <- c_convolutionEncode(message,
                              convEncoder$N,
                              convEncoder$M,
                              convEncoder$nextState,
                              convEncoder$output);

  return(code);
}

conv_decode <- function(code, convEncoder) {
   
   output <- c_convolutionDecode(code,
                                 convEncoder$N,
                                 convEncoder$M,
                                 convEncoder$prevState,
                                 convEncoder$output);
   
   # output = list(softOutput, hardOutput)
   return(output);
}

conv_decode_hard <- function(code, convEncoder) {

   output <- c_convolutionDecode_hard(code,
                                      convEncoder$N,
                                      convEncoder$M,
                                      convEncoder$prevState,
                                      convEncoder$output);

   return(output);
}