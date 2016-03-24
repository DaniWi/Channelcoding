dyn.load(paste0("src/convolution", .Platform$dynlib.ext))

generateMatrices <- function(N, M, generators) {
   
   # vector lengths of returned matrices
   num_states <- 2^M;
   matrix_dim <- 2 * num_states;
   previous_dim <- 2 * 2 * num_states;
   
   result <- .C("wrapper_generateMatrices",
                as.integer(N),
                as.integer(M),
                as.integer(generators),
                nextState = integer(length=matrix_dim),
                prevState = integer(length=previous_dim),
                output = integer(length=matrix_dim));
   
   convEncoder <- list(N = N,
                       M = M,
                       generators = generators,
                       nextState = result$nextState,
                       prevState = result$prevState,
                       output = result$output);
   
   return(convEncoder);
}

convEncode <- function(input, convEncoder) {
   # length of input message
   inputLength <- length(input);
   
   # length of resulting code
   # input bits + M termination bits, N output bits per input bit
   codeLength <- (inputLength + convEncoder$M) * convEncoder$N;
   
   result <- .C("wrapper_convolution_encode",
                as.integer(input),
                as.integer(inputLength),
                as.integer(convEncoder$N),
                as.integer(convEncoder$M),
                as.integer(convEncoder$nextState),
                as.integer(convEncoder$output),
                code = integer(length=codeLength));
   
   return(result$code);
}

convDecode <- function(code, convEncoder) {
   
   codeLength <- length(code);
   num_states <- 2^convEncoder$M;
   outputLength <- codeLength / convEncoder$N;
   
   result <- .C("wrapper_convolution_decode",
                as.integer(code),
                as.integer(codeLength),
                as.integer(convEncoder$N),
                as.integer(convEncoder$M),
                as.integer(num_states),
                as.integer(convEncoder$nextState),
                as.integer(convEncoder$prevState),
                as.integer(convEncoder$output),
                softOutput = double(outputLength),
                hardOutput = integer(outputLength));
   
   mask <- c(1:(outputLength - convEncoder$M));
   return( list(soft=result$softOutput[mask], hard=result$hardOutput[mask]));
}