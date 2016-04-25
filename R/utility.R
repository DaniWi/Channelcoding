#' @export
GetPuncturingMatrix <- function(puncturing.vector, coder.info) {

  if (any((puncturing.vector != 1)[puncturing.vector != 0])) {
    # puncturing.vector has elements with value different from 0 or 1 which is not allowed
    stop("Invalid puncturing vector! Only values 0/1 are allowed!")
  }

  if (is.null(coder.info$N)) {
    stop("Encoder has to specify list element N!")
  }

  mat <- matrix(puncturing.vector, nrow = coder.info$N)

  if (any(colSums(mat) == 0)) {
    stop("Puncturing Matrix has a column only with 0, this is illegal!")
  }

  return(mat)
}

#' @export
PunctureCode <- function(original.code, puncturing.matrix) {
  mask <- as.logical(puncturing.matrix)

  if(length(original.code) < length(mask)) {
    return(original.code[mask[1:length(original.code)]])
  } else {
    return(original.code[mask])
  }
}

#' @export
InsertPuncturingBits <- function(punctured.code, puncturing.matrix) {
  result <- c_insert_puncturing_bits(punctured.code, as.numeric(puncturing.matrix), (dim(puncturing.matrix))[1], (dim(puncturing.matrix))[2])
  return(result)
}

isOctal <- function(generators) {
  # regex check for octal number format
  x <- regexpr("^[1-7][0-7]*$",generators)

  if (length(x[x < 0]) > 0) {
    # there are numbers that are NOT in octal form
    return(FALSE)
  }

  return(TRUE)
}

maskGenerators <- function(generators, max.generator.octal) {

  new.generators = bitwAnd(max.generator.octal, generators)

  return(new.generators)
}

isCatastrophicEncoder <- function(generators) {

  # convert octal generators to decimal
  generators <- sapply(generators,octalToDecimal)

  # get GCD of all generators
  gcd <- numbers::mGCD(generators)

  # if gcd is a power of 2 the encoder is not catastrophic
  # this is done by bit checking
  return(!bitwAnd(gcd, gcd - 1) == 0)
}

octalToDecimal <- function(x) {
  dec <- 0
  i <- 0
  while (x > 0) {
    dec <- dec + (x %% 10)*(8^i)
    i <- i + 1
    x = x %/% 10
  }

  return(dec)
}

decimalToOctal <- function(x) {
  oct <- 0
  i <- 0
  while (x > 0) {
    oct <- oct + (x %% 8)*(10^i)
    i <- i + 1
    x = x%/% 8
  }

  return(oct)
}

Interleave <- function(v1, v2, v3) {
  ord1 <- 3*(1:length(v1)) - 2
  ord2 <- 3*(1:length(v2)) - 1
  ord3 <- 3*(1:length(v3))
  return(c(v1,v2,v3)[order(c(ord1,ord2,ord3))])
}

Deinterleave <- function(vec, index) {
  index <- c(rep(FALSE,index - 1), TRUE, rep(FALSE, 3 - index))
  return(vec[index])
}
