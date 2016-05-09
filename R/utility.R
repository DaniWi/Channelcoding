#' get punctuation matrix from vector
#'
#' creates a punctuation matrix from the passed punctuation vector and the
#' passed coder
#' @param punctuation.vector vector containing the punctuation information which will
#'     be transformed to a punctuation matrix
#' @param coder.info channelcoder which is used for the matrix dimensions
#' @return punctuation matrix suitable for encode and decode
#' @export
GetPunctuationMatrix <- function(punctuation.vector, coder.info) {

  if (any((punctuation.vector != 1)[punctuation.vector != 0])) {
    # punctuation.vector has elements with value different from 0 or 1 which is not allowed
    stop("Invalid punctuation vector! Only values 0/1 are allowed!")
  }

  if (is.null(coder.info$N)) {
    stop("Encoder has to specify list element N!")
  }

  if (length(punctuation.vector) %% coder.info$N != 0) {
    stop("Wrong length of punctuation vector! Must be a multiple of N (amount of exits)!")
  }

  mat <- matrix(punctuation.vector, nrow = coder.info$N)

  if (any(colSums(mat) == 0)) {
    stop("Punctuation matrix should not have a 0 column!")
  }

  return(mat)
}

PunctureCode <- function(original.code, punctuation.matrix) {
  mask <- as.logical(punctuation.matrix)

  if(length(original.code) < length(mask)) {
    punctured.code <- original.code[head(mask, length(original.code))]
    return(punctured.code)
  } else {
    return(original.code[mask])
  }

}

InsertPunctuationBits <- function(punctured.code, punctuation.matrix) {
  rows <- nrow(punctuation.matrix)
  cols <- ncol(punctuation.matrix)
  result <- c_insert_punctuation_bits(punctured.code, as.numeric(punctuation.matrix), rows, cols)
  return(result)
}

IsOctal <- function(generators) {
  # regex check for octal number format
  x <- regexpr("^[1-7][0-7]*$", generators)

  if (length(x[x < 0]) > 0) {
    # there are numbers that are NOT in octal form
    return(FALSE)
  }

  return(TRUE)
}

IsCatastrophicEncoder <- function(generators) {
  # convert octal generators to decimal
  generators <- sapply(generators, OctalToDecimal)

  # get GCD of all generators (C++ function)
  gcd <- gcd_polynomial(generators)

  # if gcd is a power of 2 the encoder is not catastrophic
  # this is done by bit checking
  return(!bitwAnd(gcd, gcd - 1) == 0)
}

CheckCoder <- function(coder) {
  minimum.fields <-
    c("N", "M", "generators", "next.state", "prev.state", "output", "rsc", "termination")

  if (!all(minimum.fields %in% names(coder))) {
    stop("Coder does not contain all fields!")
  }

  for (field in minimum.fields) {
    if (is.null(coder[[field]])) {
      stop("NULL element in the coder!")
    }
  }

  stopifnot(coder$N > 1, coder$M > 0)

  if (length(coder$generators) < coder$N) {
    stop("Coder has to few generator polynoms!")
  }

  if (!IsOctal(coder$generators)) {
    stop("At least one of the generators is not in octal form!")
  }

  max.generator.octal = DecimalToOctal(2^(coder$M+1) - 1)

  if (any(coder$generators > max.generator.octal)) {
    stop("At least one of the generators is bigger than allowed!")
  }
}

OctalToDecimal <- function(x) {
  dec <- 0
  i <- 0
  while (x > 0) {
    dec <- dec + (x %% 10) * (8^i)
    i <- i + 1
    x = x %/% 10
  }

  return(dec)
}

DecimalToOctal <- function(x) {
  oct <- 0
  i <- 0
  while (x > 0) {
    oct <- oct + (x %% 8) * (10^i)
    i <- i + 1
    x = x %/% 8
  }

  return(oct)
}

Interleave <- function(v1, v2, v3) {
  ord1 <- 3*(1:length(v1)) - 2
  ord2 <- 3*(1:length(v2)) - 1
  ord3 <- 3*(1:length(v3))
  return(c(v1,v2,v3)[order(c(ord1, ord2, ord3))])
}

Deinterleave <- function(vec, index) {
  index <- c(rep(FALSE,index - 1), TRUE, rep(FALSE, 3 - index))
  return(vec[index])
}

Shift <- function(vector, n) {
  l <- length(vector)
  if (n %% l == 0) {
    return(vector)
  }

  result <- c(tail(vector, n %% l), head(vector, (l - n) %% l))
}

MaskGenerators <- function(generators, max.generator.octal) {

  new.generators = bitwAnd(max.generator.octal, generators)

  return(new.generators)
}


