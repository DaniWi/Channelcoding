# interface for convolutional codes
# provided functions:
#  - generate nsc coder
#  - generate rsc coder
#  - encode
#  - decode (soft in & out)
#  - decode (hard decision)

library(Rcpp);
sourceCpp('src/convolution.cpp');

# generate nsc coder
generateConvEncoder_nsc <- function(N, M, generators) {
   
   # nsc requires N generator polynoms
   if (length(generators) < N) {
      # stop execution if too few generators
      stop("Too few generator polynoms!");
   }
   else if (length(generators) > N) {
      # just use first N polynoms if too many are provided
      generators <- head(generators, N);
   }
   
   matrixList <- c_generateMatrices_nsc(N,M,generators);

   convEncoder <- list(N = N,
                       M = M,
                       generators = generators,
                       nextState = matrixList$nextState,
                       prevState = matrixList$prevState,
                       output = matrixList$output,
                       nsc = TRUE,
                       termination = NULL);

   return(convEncoder);
}

# generate rsc encoder
generateConvEncoder_rsc <- function(N, M, generators) {
   
   # rsc requires N+1 generator polynoms
   if (length(generators) < N+1) {
      # stop execution if too few generators
      stop("Too few generator polynoms!");
   }
   else if (length(generators) > N+1) {
      # just use first N polynoms if too many are provided
      generators <- head(generators, N+1);
   }
   
   matrixList <- c_generateMatrices_rsc(N,M,generators);
   
   convEncoder_rsc <- list(N = N,
                           M = M,
                           generators = generators,
                           nextState = matrixList$nextState,
                           prevState = matrixList$prevState,
                           output = matrixList$output,
                           nsc = FALSE,
                           termination = matrixList$termination);
   
   return(convEncoder_rsc);
}

# encode of message with convolutional encoder
conv_encode <- function(message, convEncoder, terminate = TRUE) {

   code <- c_convolutionEncode(message,
                               convEncoder$N,
                               convEncoder$M,
                               convEncoder$nextState,
                               convEncoder$output,
                               as.integer(convEncoder$nsc),
                               convEncoder$termination,
                               as.integer(terminate));
   return(code);
}

conv_decode <- function(code, convEncoder, terminate = TRUE) {
   
   output <- c_convolutionDecode(code,
                                 convEncoder$N,
                                 convEncoder$M,
                                 convEncoder$prevState,
                                 convEncoder$output);
   
   # if terminated, termination bits are thrown away
   if (terminate == TRUE) {
      soft <- head(output$softOutput, length(output$softOutput) - N*M);
      hard <- head(output$hardOutput, length(output$hardOutput) - N*M);
      
      newlist <- list(softOutput = soft, hardOutput = hard);
      return(newlist);
   }
   
   return(output);
}

conv_decode_hard <- function(code, convEncoder, terminate = TRUE) {

   output <- c_convolutionDecode_hard(code,
                                      convEncoder$N,
                                      convEncoder$M,
                                      convEncoder$prevState,
                                      convEncoder$output);

   # if terminated, termination bits are thrown away
   if (terminate == TRUE) {
      return(head(output, length(output) - N*M));
   }
   
   return(output);
}