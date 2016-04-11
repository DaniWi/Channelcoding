isOctal <- function(generators) {
  # regex check for octal number format
  x <- regexpr("^[1-7][0-7]*$",generators)

  if (length(x[x < 0]) > 0) {
    # there are numbers that are NOT in octal form
    return(FALSE)
  }

  return(TRUE)
}

maskGenerators <- function(generators, M) {
  max.generator = 2^M - 1

  new.generators = bitwAnd(max.generator, generators)

  return(new.generators)
}


isCatastrophicEncoder <- function(generators) {

  return(FALSE)
}
