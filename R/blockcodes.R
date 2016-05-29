# blockcodes.R
#
# interface for block codes
# provided functions:
#  - generate hamming encoder
#  - generate bch encoder
#  - encode
#  - decode
#  - simulation
#  - open pdf



#' Generate Hamming encoder.
#'
#' Generates a block encoder for Hamming block codes.
#' @details The resulting Hamming Code is obtained by constructing a (code.length - data.length) x code.length Parity-Check-Matrix.
#' The columns are the binary vectors of the integers 1 to code.length. Column Permutations are then applied to get the
#' Form (t(A) | Ir). The resulting Generator Matrix is (In | A).
#' @param code.length Length of the encoded message
#' @param data.length Length of the plain message
#' @return A block encoder represented as a list containing:
#'     type, code.length, data.length,
#'     2 Matrices gen.matrix and check.matrix, being the Generator Matrix and the Parity Check Matrix
#' @examples
#' # standard Hamming encoder with code-rate = 0.57
#' BlockGenerateEncoderHamming(7,4)
#' @author Benedikt Wimmer
#' @export
BlockGenerateEncoderHamming = function(code.length = 7, data.length = 4){

  binVector <- function(n, min.length) {
    result <- integer()
    i <- 0
    while (n != 0) {
      result <- c(n %% 2, result)
      n = n %/% 2
      i <- i + 1
    }
    zeros <- rep(0, times = min.length)
    result <- c(zeros, result)
    return(tail(result,min.length))
  }

  r = code.length - data.length

  if((2**r-1) != code.length) {
    stop("Bad Parameters! '2^(code.length - data.length) - 1 = code.length' has to be satisfied, eg (7,4), (15,11),etc.")
  }

  p = 1:code.length
  p = p[(log2(p)%%1)!=0 ]
  p = sapply(p, function(x) binVector(x,r))
  p = matrix(p, nrow = r)

  check.matrix = cbind(p,diag(r))
  gen.matrix = cbind(diag(data.length),t(p))

  block.encoder = list(type = 'HAMMING',
                       code.length = code.length,
                       data.length = data.length,
                       gen.matrix = gen.matrix,
                       check.matrix = check.matrix
  )
  return(block.encoder)
}


#' Generate BCH encoder.
#'
#' Generates a block encoder for BCH block codes.
#' @details Apart from necessary parameters like data.length, a BCH encoder consists of a Generator Polynomial, which is needed for encoding,
#' and the polynomial and index representation of the underlying Galois Field GF(code.length+1), which is needed for constructing and solving the
#' error location Polynomial.
#' @param code.length Length of the encoded message
#' @param code.t Error correcting capability
#' @return A block encoder represented as a list containing:
#'     type, code.length, code.m(minimal m for 2^m >= code.length), data.length, code.t,
#'     gen.poly,
#'     2 Vectors alpha_to and index_of describing the Galois Field GF(2^m)
#' @examples
#' # standard bch encoder with code-rate = 0.33
#' BlockGenerateEncoderBCH(15,3)
#' @author Benedikt Wimmer
#' @export
BlockGenerateEncoderBCH = function(code.length = 15, code.t = 3){

  if(2*code.t > code.length)
    stop("Invalid Parameters for BCH, 2*t must be < code.length")

  m = ceiling(log2(code.length))
  ret = c_getGeneratorPoly(code.length,
                           m,
                           code.t)

  block.encoder = list(type = 'BCH',
                       code.length = code.length,
                       code.m = m,
                       data.length = ret$data.length,
                       code.t = code.t,
                       gen.poly = ret$gen.poly,
                       alpha_to = ret$alpha_to,
                       index_of = ret$index_of
                       )
  return(block.encoder)

}

#' Block encoding of a message.
#'
#' Produces a block code of a message.
#' @param message The message to be encoded.
#' @param block.encoder Block encoder used for encoding.
#' @param visualize If TRUE a beamer PDF file is generated showing the encode process.
#' @return The encoded message as binary vector.
#' @examples
#' coder <- BlockGenerateEncoderBCH(15,3)
#' plain <- c(1,0,0,1,1)
#'
#' # standard encoding
#' BlockEncode(plain, coder)
#'
#' # standard encoding with visualization
#' # BlockEncode(plain, coder, visualize = TRUE)
#'
#' # use default values
#' BlockEncode(plain)
#' @author Benedikt Wimmer
#' @export
BlockEncode = function(message, block.encoder= NULL, visualize=FALSE){

  stopifnot(length(message) > 0)

  if (any((message != 1)[message != 0])) {
    stop("Message should only contain 0 and 1!")
  }

  if (is.null(block.encoder)) {
    warning("Standard Block Encoder was used! BCH(15,5,7)")
    block.encoder <- BlockGenerateEncoderBCH()
  }

  if((length(message)%%block.encoder$data.length)!= 0){
    message = c(message, rep(0,(block.encoder$data.length-(length(message)%%block.encoder$data.length))));
  }

  msgMatrix = t(matrix(message, nrow=block.encoder$data.length));

  switch(block.encoder$type,
         HAMMING={
           retMatrix = apply(msgMatrix,1,function(x) MatrixEncode(x, block.encoder))
           if(visualize){

             if (block.encoder$code.length > 25) {
               warning("Code.length is too long, matrices can't be displayed properly. No PDF was created.")
             }else{
             rmarkdown::render(system.file("rmd", "BlockEncodeHamming.Rmd", package = "channelcoding"),
                               output_dir = system.file("pdf", package = "channelcoding"),
                               encoding = "UTF-8",
                               params = list(block.encoder = block.encoder,
                                             message = message,
                                             code = as.vector(retMatrix)))

             rstudioapi::viewer(system.file("pdf", "BlockEncodeHamming.pdf", package = "channelcoding"))
             }
           }
         },
         BCH={



           retMatrix = apply(msgMatrix,1,function(x) c_bchEncode(x, block.encoder$gen.poly, block.encoder$code.length, block.encoder$data.length))
           if(visualize){

             if (block.encoder$code.length > 25) {
               warning("Code.length is too long, polynomials can't be displayed properly. No PDF was created.")
             }else{
             rmarkdown::render(system.file("rmd", "BlockEncodeBCH.Rmd", package = "channelcoding"),
                               output_dir = system.file("pdf", package = "channelcoding"),
                               encoding = "UTF-8",
                               params = list(block.encoder = block.encoder,
                                             message = message,
                                             code = as.vector(retMatrix)))

             rstudioapi::viewer(system.file("pdf", "BlockEncodeBCH.pdf", package = "channelcoding"))
             }
           }
           },
         {
           print("No or Wrong Type selected")
            return()
         })

  return(as.vector(retMatrix))
}


#' Block decoding of a code.
#'
#' Decodes a block codeword.
#' @param code The code to be decoded.
#' @param block.encoder Block encoder used for encoding.
#' @param visualize If TRUE a beamer PDF file is generated showing the decode process.
#' @return The decoded message as binary Vector
#' @examples
#' coder <- BlockGenerateEncoderBCH(15,3)
#' plain <- c(1,0,0,1,1)
#'
#' # standard encoding and decoding
#' coded <- BlockEncode(plain, coder)
#' BlockDecode(coded, coder)
#'
#' # with visualization
#' # coded <- BlockEncode(plain, coder)
#' # BlockDecode(coded, coder, visualize = TRUE)
#'
#' # with message distortion
#' coded <- BlockEncode(plain, coder)
#' noisy <- ApplyNoise(coded, SNR.db = 3, binary = TRUE)
#' BlockDecode(noisy, coder)
#' @author Benedikt Wimmer
#' @export
BlockDecode = function(code, block.encoder = NULL, visualize=FALSE){

  stopifnot(length(code) > 0)

  if (any((code != 1)[code != 0])) {
    stop("Message should only contain 0 and 1!")
  }

  if (is.null(block.encoder)) {
    warning("Standard Block Encoder was used! BCH(15,5,7)")
    block.encoder <- BlockGenerateEncoderBCH()
  }

  if((length(code)%%block.encoder$code.length)!= 0){
    code = c(code, rep(0,(block.encoder$code.length-(length(code)%%block.encoder$code.length))));
  }

  codeMatrix = t(matrix(code, nrow=block.encoder$code.length));

  switch(block.encoder$type,
         HAMMING={
           retListMatrix = apply(codeMatrix,1,function(x) MatrixDecode(x, block.encoder))
           retMatrix = sapply(retListMatrix,function(x) x$decoded)

            if(visualize){
              if (block.encoder$code.length > 25) {
                warning("Code.length is too long, matrices can't be displayed properly. No PDF was created.")
              }else{
               correctedMatrix = sapply(retListMatrix,function(x) x$input.corrected)

               rmarkdown::render(system.file("rmd", "BlockDecodeHamming.Rmd", package = "channelcoding"),
                                 output_dir = system.file("pdf", package = "channelcoding"),
                                 encoding = "UTF-8",
                                 params = list(block.encoder = block.encoder,
                                               enc = code,
                                               enc.corrected = as.vector(correctedMatrix),
                                               dec = as.vector(retMatrix)))

               rstudioapi::viewer(system.file("pdf", "BlockDecodeHamming.pdf", package = "channelcoding"))
              }
            }
         },
         BCH={
           retListMatrix = apply(codeMatrix,1,function(x) c_bchDecode(x,
                                                                 block.encoder$code.length,
                                                                 block.encoder$code.m,
                                                                 block.encoder$data.length,
                                                                 block.encoder$code.t,
                                                                 block.encoder$alpha_to,
                                                                 block.encoder$index_of))

           retMatrix = sapply(retListMatrix,function(x) x$decoded)

            if(visualize){
              if (block.encoder$code.length > 25) {
                warning("Code.length is too long, polynomials can't be displayed properly. No PDF was created.")
              }else{

              correctedMatrix = sapply(retListMatrix,function(x) x$input.corrected)
              failedMatrix = sapply(retListMatrix,function(x) x$failed)
              rmarkdown::render(system.file("rmd", "BlockDecodeBCH.Rmd", package = "channelcoding"),
                                output_dir = system.file("pdf", package = "channelcoding"),
                                encoding = "UTF-8",
                                params = list(block.encoder = block.encoder,
                                              enc = code,
                                              enc.corrected = as.vector(correctedMatrix),
                                              dec = as.vector(retMatrix),
                                              failed = as.vector(failedMatrix)))

              rstudioapi::viewer(system.file("pdf", "BlockDecodeBCH.pdf", package = "channelcoding"))
              }
            }

           },
         {
           print("No or Wrong Type selected")
            return()
         })

  return(as.vector(retMatrix))

}

#' Block Simulation.
#'
#' Simulation of a block encode and decode process over a noisy channel.
#'
#' @param coder Block coder used for the simulation. Can be created via
#'     \code{\link{BlockGenerateEncoderBCH}} or \code{\link{BlockGenerateEncoderHamming}}.
#' @param msg.length Message length of the randomly created messages to be encoded. Will be shortened to a multiple of the coders blocksize
#' @param min.db Minimum SNR to be tested.
#' @param max.db Maximum SNR to be tested.
#' @param db.interval Step between two SNRs tested.
#' @param iterations.per.db Number of encode and decode processes per SNR.
#' @param visualize If true a PDF report is generated.
#' @return Dataframe containing the bit-error-rates for each SNR tested.
#' @examples
#' # use all default parameters
#' BlockSimulation()
#'
#' # Custom coder
#' coder <- BlockGenerateEncoderHamming(15,11)
#' BlockSimulation(coder, 15, 0.01, 1, 0.05, 50, FALSE)
#' @author Benedikt Wimmer
#' @export
BlockSimulation <- function(coder = NULL,
                           msg.length = 100,
                           min.db = 0.1,
                           max.db = 2.0,
                           db.interval = 0.1,
                           iterations.per.db = 100,
                           visualize = FALSE)
{
  stopifnot(msg.length > 0, iterations.per.db > 0,
            min.db > 0, max.db > 0, max.db >= min.db, db.interval > 0)

  if (is.null(coder)) {
    warning("Standard Block Encoder was used! BCH(15,5,7)")
    coder <- BlockGenerateEncoderBCH(15,3)
  }

  #msg.length an Blockgröße anpassen, sonst haben decoded und message evt. untersch. Längen
  if(msg.length%%coder$data.length != 0){
    msg.length = msg.length - (msg.length%%coder$data.length)
  }


  v.db <- seq(from = min.db, to = max.db, by = db.interval)
  v.ber <- numeric(0)

  total.errors <- 0

  total.iterations <- length(v.db) * iterations.per.db

  print("Block Simulation")
  progress.bar <- txtProgressBar(min = 0, max = total.iterations, style = 3)

  progress.counter <- 0

  for (db in v.db) {
    for (i in 1 : iterations.per.db) {
      # message erzeugen
      message <- sample(c(0,1), msg.length, replace = TRUE)

      # encode
      coded <- BlockEncode(message, coder)

      # noise hinzufügen
      noisy <- ApplyNoise(coded, db, binary = TRUE)

      # anzahl flipped bits (channel errors)
      # coded.hard <- ifelse(coded >= 0, 0, 1)
      # noisy.hard <- ifelse(noisy >= 0, 0, 1)
      # channel.errors <- sum(abs(coded.hard - noisy.hard))

      # decode
      decoded <- BlockDecode(noisy, coder)

      # vgl decoded & message
      decode.errors <- sum(abs(decoded - message))

      total.errors <- total.errors + decode.errors

      progress.counter <- progress.counter + 1

      setTxtProgressBar(progress.bar, progress.counter)
    }

    v.ber <- c(v.ber, total.errors / (msg.length * iterations.per.db))
    total.errors <- 0
  }

  close(progress.bar)

  df <- data.frame(db = v.db, ber = v.ber)

  if (visualize) {
    rmarkdown::render(system.file("rmd", "BlockSimulation.Rmd", package = "channelcoding"),
                      output_dir = system.file("pdf", package = "channelcoding"),
                      output_file = "BlockSimulation.pdf",
                      encoding = "UTF-8",
                      params = list(message.length = msg.length,
                                    iterations.per.db = iterations.per.db,
                                    min.db = min.db,
                                    max.db = max.db,
                                    db.interval = db.interval,
                                    encoder = coder,
                                    dataframe = df))

    rstudioapi::viewer(system.file("pdf", "BlockSimulation.pdf", package = "channelcoding"))
  }

  return(df)
}


#' Open visualization PDF.
#'
#' With this function it is easy to reopen the PDF files created with
#'     \code{\link{BlockEncode}}, \code{\link{BlockDecode}},
#'     and \code{\link{BlockSimulation}}.
#'     The files are stored in the program files of R.
#'     If the corresponding file does not exist yet an error message is printed.
#' @param encode Flag to open encode PDFs (if true) or decode PDFs (if false).
#' @param hamming Flag to open encode or decode PDFs for Hamming codes. For BCH use false.
#' @param simulation Flag to open simulation PDFs. This flag has highest precedence
#'     meaning that if the simulation flag is TRUE the function will look for the
#'     simulation PDF. Only if FALSE the others are evaluated.
#' @examples
#' # open encode for Hamming PDF
#' BlockOpenPDF(encode = TRUE, hamming = TRUE)
#'
#' # open decode for BCH PDF
#' BlockOpenPDF(encode = FALSE, hamming = FALSE)
#'
#' # open block simulation PDF
#' BlockOpenPDF(simulation = TRUE)
#' @author Benedikt Wimmer
#' @export
BlockOpenPDF <- function(encode = TRUE, hamming = FALSE, simulation = FALSE) {
  if (simulation) {
    path <- system.file("pdf", "BlockSimulation.pdf", package = "channelcoding")
  } else {
    if (encode) {
      if (hamming) {
        path <- system.file("pdf", "BlockEncodeHamming.pdf", package = "channelcoding")
      } else {
        path <- system.file("pdf", "BlockEncodeBCH.pdf", package = "channelcoding")
      }
    } else {
      if (hamming) {
        path <- system.file("pdf", "BlockDecodeHamming.pdf", package = "channelcoding")
      } else {
        path <- system.file("pdf", "BlockDecodeBCH.pdf", package = "channelcoding")
      }
    }
  }
  if (path != "") {
    rstudioapi::viewer(path)
  } else {
    warning("File does not exists!")
  }
}



MatrixEncode = function(message, block.encoder)
{
  return((message%*%block.encoder$gen.matrix)%%2)
}

MatrixDecode = function(message, block.encoder)
{
  result = as.vector((block.encoder$check.matrix %*% message) %% 2)

  FindErrors = function(){
     i = 1
    function(x){
       if(isTRUE(identical(result,x))){
              message[i] <<- (message[i]+1)%%2
       }
      i <<- i+1
      x
    }
  }

apply(block.encoder$check.matrix, 2, FindErrors())

  ret = list(decoded = message[1:block.encoder$data.length],
             input.corrected = message)

  return(ret)
}

