#' Message distortion.
#'
#' This function will randomly distort a clean message.
#' @author Bene Wimmer
#' @param msg Message to alter.
#' @param SNR.db Signal-Noise-Ratio in dB simulating the noisy channel.
#' @param binary False = Soft Output, True = Hard Output
#'
#' @return Distorted message containing noise.
#'
#' @export
ApplyNoise <- function(msg, SNR.db = 3, binary = FALSE)
{
  msg.len <- length(msg);
  SNR.linear <- 10^(SNR.db/10);
  power <- sum(msg^2)/(msg.len); #power of vector msg

  # noise is a vector of lenth msg_len and contains
  # normal distributed values with mean 0 and standard-deviation 1
  noise <- sqrt(power / SNR.linear) * rnorm(msg.len,0,1)

  msg.out <- msg + noise

  if(binary){
    msg.out <- ifelse(msg.out <= 0.5, 0, 1)
  }

  return(msg.out);
}

#' Channelcoding Simulation
#'
#' Simulation of channelcoding techniques (blockcodes, convolutional codes
#'     and turbo codes) and comparison of their bit-error-rates.
#'
#' @param msg.length Message length of the randomly created messages to be encoded.
#' @param min.db Minimum SNR to be tested.
#' @param max.db Maximum SNR to be tested.
#' @param db.interval Step between two SNRs tested.
#' @param iterations.per.db Amount of iterations each signal/noise ration step.
#' @param turbo.decode.iterations Amount of decoding iterations inside the turbo decoder.
#' @param visualize If true a PDF report is generated.
#'
#' @return Dataframe containing the bit-error-rates for each coding technique
#'     and each SNR tested.
#'
#' @export
ChannelcodingSimulation <- function(msg.length = 100,
                                    min.db = 0.1,
                                    max.db = 2.0,
                                    db.interval = 0.1,
                                    iterations.per.db = 100,
                                    turbo.decode.iterations = 5,
                                    visualize = FALSE)
{
  block.df <- BlockSimulation(msg.length = msg.length,
                              min.db = min.db,
                              max.db = max.db,
                              db.interval = db.interval,
                              iterations.per.db = iterations.per.db)

  conv.df <- ConvSimulation(msg.length = msg.length,
                            min.db = min.db,
                            max.db = max.db,
                            db.interval = db.interval,
                            iterations.per.db = iterations.per.db)

  turbo.df <- TurboSimulation(msg.length = msg.length,
                              min.db = min.db,
                              max.db = max.db,
                              db.interval = db.interval,
                              iterations.per.db = iterations.per.db,
                              decode.iterations = turbo.decode.iterations)

  df <- data.frame(db = block.df$db,
                   block.ber = block.df$ber,
                   conv.ber = conv.df$ber,
                   turbo.ber = turbo.df$ber)


  if (visualize) {
    rmarkdown::render(system.file("rmd", "ChannelcodingSimulation.Rmd", package = "channelcoding"),
                      output_dir = system.file("pdf", package = "channelcoding"),
                      output_file = "ChannelcodingSimulation.pdf",
                      encoding = "UTF-8",
                      params = list(message.length = msg.length,
                                    min.db = min.db,
                                    max.db = max.db,
                                    db.interval = db.interval,
                                    dataframe = df))

    rstudioapi::viewer(system.file("pdf", "ChannelcodingSimulation.pdf", package = "channelcoding"))
  }

  return(df)
}

#' Show simulation data in a plot
#'
#' Shows passed simulation data (can be created by the functions
#' \code{\link{ConvSimulation}}, \code{\link{TurboSimulation}}
#' and \code{\link{BlockSimulation}}) containing the bit-error-rates for several SNRs
#' in one plot for easy comparison.
#' @param ... Dataframes created by the simulation functions.
#' @examples
#' # create dataframes from simulation
#' block <- BlockSimulation()
#' conv <- ConvSimulation()
#' turbo <- TurboSimulation()
#'
#' # show in plot
#' PlotSimulationData(block, conv, turbo)
#' @export
PlotSimulationData <- function(...) {
  arguments <- list(...)
  if (!all(sapply(arguments, is.data.frame))) {
    stop("Arguments are not data frames!")
  }
  if (any(sapply(arguments, function(x) is.null(x[["ber"]]) || is.null(x[["db"]])))) {
    stop("The data frame are not correct! (required columns: ber, db)")
  }
  min.ber <- min(sapply(arguments, function(x) min(x[["ber"]])))
  max.ber <- max(sapply(arguments, function(x) max(x[["ber"]])))

  new.df <- Reduce(rbind, arguments)
  new.df$Arguments <- rep(paste("Argument", 1:length(arguments)), sapply(arguments, nrow))
  ggplot2::ggplot(new.df, ggplot2::aes(new.df[["db"]], new.df[["ber"]], color = new.df[["Arguments"]])) +
    ggplot2::geom_line() +
    ggplot2::xlab("Signal Rausch Verh\u00e4ltnis [dB]") +
    ggplot2::ylab("Bitfehlerrate") +
    ggplot2::theme(legend.title = ggplot2::element_blank())
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

IsCatastrophicEncoder <- function(generators.oct, M) {
  # convert octal generators to decimal
  generators.dec <- sapply(generators.oct, OctalToDecimal)

  # turn bits of generators
  generators.dec <- sapply(generators.dec, TurnBitsRound, nbits = M)

  # get GCD of all generators (C++ function)
  gcd <- gcd_polynomial(generators.dec)

  # if gcd is a power of 2 the encoder is not catastrophic
  # this is done by bit checking
  return(!bitwAnd(gcd, gcd - 1) == 0)
}

TurnBitsRound <- function(number, nbits) {
  result <- 0

  for (i in 0:(nbits-1)) {
    bit <- bitwAnd(bitwShiftR(number, i), 1)
    bit <- bitwShiftL(bit, nbits - i - 1)
    result <- bitwOr(result, bit)
  }

  return(result)
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


