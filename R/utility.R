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
