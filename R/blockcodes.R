.HammingEncode = function(msg, params, visualize)
{
  if(length(params)!= 2){
    print("Wrong Parameters for type HAMMING")
    return()
  }
  r = params[1]-params[2];
  p = apply(apply(diag(r),1,function(x) x+1),1, function(x) x%%2);


  p = rbind(p, matrix(rep(1,params[2]-r), ncol = r));

  genMatrix = cbind(diag(params[2]), p);

  if(visualize){
    print("Generator Matrix:")
    print(genMatrix)
  }

  if((length(msg)%%params[2])!= 0){
    msg = c(msg, rep(0,(params[2]-(length(msg)%%params[2]))));
  }
  msgMatrix = t(matrix(msg, nrow=params[2]));
  retMatrix = apply(msgMatrix,1,function(x) .MatrixEncode(x, genMatrix, visualize))

  return(as.vector(retMatrix))

}

.HammingDecode = function(msg, params, visualize)
{

  if(length(params)!= 2){
    print("Wrong Parameters for type HAMMING")
    return()
  }

  r = params[1]-params[2];

  p = apply(apply(diag(r),1,function(x) x+1),1, function(x) x%%2);
  p = rbind(p, matrix(rep(1,params[2]-r), ncol = r));

  checkMatrix = cbind(t(p),diag(r));


  if((length(msg)%%params[1])!= 0){
    msg = c(msg, rep(0,(params[1]-(length(msg)%%params[1]))));
  }

  decmsg = c();
  for(k in 1:(length(msg)/params[1])){

    subMsg = msg[((k-1)*params[1]+1):((k-1)*params[1]+params[1])];
    result = .MatrixDecode(subMsg, checkMatrix, visualize)[,1];

    errorIndex = -1;
    for(i in (1:params[1])){
      if(isTRUE(all.equal(result,checkMatrix[,i]))){
        errorIndex = i;
      }
    }


    if(errorIndex>0){
      subMsg[errorIndex] = (subMsg[errorIndex+1])%%2;
    }

    decmsg = c(decmsg,subMsg[1:params[2]])
  }

  if(visualize){
    print("Checkmatrix:")
    print(checkMatrix);
  }

  return(decmsg)

}


.MatrixEncode = function(msg, genMatrix, visualize)
{

  if(length(msg)!= nrow(genMatrix)){
    print("Wrong dimensions!")
    return()
  }

  return((msg%*%genMatrix)%%2)

}

.MatrixDecode = function(msg, checkmatrix, visualize)
{
  if(length(msg)!= ncol(checkmatrix)){
    print("Wrong dimensions!")
    return()
  }

  return((checkmatrix%*%msg)%%2)
}


.BCHEncode = function(msg,params,visualize=FALSE){

  if(length(params)>2)
    print("wrong params for bch, must be c(length,t)")

  length=params[1]
  t = params[2]
  m= ceiling(log2(length))


  result = c_getGeneratorPoly(length,m,t)
  genpoly = result$genPoly
  k = result$k

  if(visualize){
    print(sprintf("This is a (%d, %d, %d) binary BCH code", length, k, 2*t+1))
    print(sprintf("Generator Polynomial g(x) = %s",paste(genpoly, collapse = "")))

  }
  # print(k)


  if((length(msg)%%k)!= 0){
    msg = c(msg, rep(0,(k-(length(msg)%%k))));
  }

  msgMatrix = t(matrix(msg, nrow=k));

  if(visualize){
    print("Computing the Codewort c(x) = c( b(x),d(x) ) with d(x) being the original message and b(x) the remainder of x^(length-k)*d(x) divided by g(x).")
    print(sprintf("Message is divided into %d block(s) of length %d",(length(msg)/k),k))
    print(apply(msgMatrix,1,function(x) sprintf("| %s |",paste(x,collapse = ""))))
  }

  retMatrix = apply(msgMatrix,1,function(x) c_bchEncode(x, genpoly, length, k))

  return(as.vector(retMatrix))

  #return(.encodePolynomial(msg,genpoly))

}

.BCHDecode = function(msg,params,visualize=FALSE){

  if(length(params)>2)
    print("wrong params for bch, must be c(length,t)")

  length=params[1]
  t = params[2]
  m= ceiling(log2(length))


  if((length(msg)%%length)!= 0){
    msg = c(msg, rep(0,(length-(length(msg)%%length))));
  }

  result = c_getGeneratorPoly(length,m,t)
  genpoly = result$genPoly
  k = result$k

  if(visualize){
    print(sprintf("This is a (%d, %d, %d) binary BCH code", length, k, 2*t+1))
    print(sprintf("Generator Polynomial g(x) = %s",paste(genpoly, collapse = "")))
  }

  #decMsg = c_bchDecode(msg,length, m, t, k, length(msg)/length)

  msgMatrix = t(matrix(msg, nrow=length));
  retMatrix = apply(msgMatrix,1,function(x) c_bchDecode(x,result$alpha_to,result$index_of,length, m, t, result$k))

  return(as.vector(retMatrix))

}

# .encodePolynomial = function(data, g, length, k, visualize=FALSE)
# {
#
# #   data_star = c(rep(0,length(g)-1),data)
# #   res= polynomial(data_star)%%polynomial(g)
# #   res = as.integer(res)
# #  res = sapply(res ,function(x) x%%2)
# #   res = as.integer(data_star - c(res,rep(0,length(data))))
# #  res = sapply(res, function(x) x%%2)
# #   return(res)
#
#   dyn.load("src/bch.so")
#   result=.C("encode_bch",length_r=as.integer(length), lm=as.integer(k), g=as.integer(g),bb=as.integer(rep(0,(length-k))), recd=as.integer(rep(0,length)), data=as.integer(data), PACKAGE = "bch")
#   if(visualize){
#     print(sprintf("x^%d*d(%s) mod g(x) = b(%s)",(length-k),paste(data, collapse = ""),paste(result$bb, collapse = "")))
#     print(sprintf("Coefficient polynomial c(x) = %s",paste(result$recd,collapse = "")))
#     }
#   return(result$recd)
#
# }

# #'@useDynLib bch
# .decodeBch = function(m,length,t,msg)
# {
#   dyn.load("src/bch.so")
#   result=.C("decode_bch",m_r=as.integer(m),length_r=as.integer(length),t_r=as.integer(t),alpha_to=as.integer(rep(0,2**m)),index_of=as.integer(rep(0,2**m)),recd=as.integer(msg), nBlocks=as.integer(length(msg)/length), PACKAGE = "bch")
#   return(result$recd)
# }
#
# #'@useDynLib bch
# .getGenPoly= function(m, length, t)
# {
#   dyn.load("src/bch.so")
#   result=.C("gen_poly",m_r=as.integer(m),length_r=as.integer(length),t_r=as.integer(t),g=as.integer(rep(0,2**m)),alpha_to=as.integer(rep(0,2**m)),index_of=as.integer(rep(0,2**m)),length_g=as.integer(0),lm=as.integer(0),PACKAGE = "bch")
#   #return(c((result$g[1:(result$length_g+1)]),(result$lm)))
#   return(result)
#   #return(result$g[1:(result$length_g+1)])
# }
#
