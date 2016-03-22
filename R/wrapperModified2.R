dyn.load(paste0("src/turbo_map_sova_modified2", .Platform$dynlib.ext))

encode <- function(input, permutation, encoder_info) {
  input_term <- c(input, rep(-1, encoder_info$amount_register))
  input_length <- length(input_term)
  output_length <- 3*input_length
  result <- .C("wrapper_encode",
               as.integer(input_term),
               out=integer(length=output_length),
               as.integer(permutation),
               as.integer(input_length),
               as.integer(output_length),
               as.integer(encoder_info$amount_register),
               as.integer(encoder_info$previous_tab),
               as.integer(encoder_info$next_tab),
               as.integer(encoder_info$parity_tab),
               as.integer(encoder_info$term_tab),
               as.integer(encoder_info$use_parity_index))
  return(result$out)
}

decode <- function(input, permutation, iterations, decoder_info){
  input_length <- length(input)
  output_length <- input_length/3
  result <- .C("wrapper_decode",
               as.double(input),
               out_soft=double(length=output_length),
               out_hard=integer(length=output_length),
               as.integer(permutation),
               as.integer(input_length),
               as.integer(output_length),
               as.integer(iterations),
               as.integer(decoder_info$amount_register),
               as.integer(decoder_info$previous_tab),
               as.integer(decoder_info$next_tab),
               as.integer(decoder_info$parity_tab),
               as.integer(decoder_info$term_tab),
               as.integer(decoder_info$use_parity_index))

  mask <- c(1:(output_length-2))
  return(list(soft = result$out_soft[mask],hard = result$out_hard[mask]))
}

get_permutation <- function(length, encoder_info){
  return(sample(c(0:(length+encoder_info$amount_register-1))))
}
