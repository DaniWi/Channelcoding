#dyn.load(paste0("src/turbo_map_sova_modified", .Platform$dynlib.ext))

encode <- function(input, permutation) {
  input_term <- c(input, -1, -1)
  input_length <- length(input_term)
  output_length <- 3*input_length
  result <- .C("wrapper_encode",
               as.integer(input_term),
               out=integer(length=output_length),
               as.integer(permutation),
               as.integer(input_length),
               as.integer(output_length))
  return(result$out)
}

decode <- function(input, permutation, iterations){
  input_length <- length(input)
  output_length <- input_length/3
  result <- .C("wrapper_decode",
               as.double(input),
               out_soft=double(length=output_length),
               out_hard=integer(length=output_length),
               as.integer(permutation),
               as.integer(input_length),
               as.integer(output_length),
               as.integer(iterations))
  mask <- c(1:(output_length-2))
  return(list(soft = result$out_soft[mask],hard = result$out_hard[mask]))
}

get_permutation <- function(length){
  return(sample(c(0:(length+1))))
}
