dyn.load(paste0("src/convolution", .Platform$dynlib.ext))

generateMatrices <- function(N, M, generators) {
   
   # vector lengths of returned matrices
   num_states <- 2^M;
   matrix_dim <- 2 * num_states;
   previous_dim <- 2 * 2 * num_states;
   
   result <- .C("generateMatrices",
                as.integer(N),
                as.integer(M),
                as.integer(generators),
                nextState = integer(length=matrix_dim),
                prevState = integer(length=previous_dim),
                output = integer(length=matrix_dim));
   
   return(result);
}

# convEncode <- function