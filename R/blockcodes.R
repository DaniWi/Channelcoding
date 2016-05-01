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


.BCHEncode = function(msg,params,visualize){

  if(length(params)>2)
    print("wrong params for bch, must be c(length,t)")

  length=params[1]
  t = params[2]
  m= ceiling(log2(length))


  result = .getGenPoly(m,length,t)
  genpoly = result$g[1:(result$length_g+1)]
  k = result$lm

 # print(k)

 if((length(msg)%%k)!= 0){
    msg = c(msg, rep(0,(k-(length(msg)%%k))));
  }

  msgMatrix = t(matrix(msg, nrow=k));
  retMatrix = apply(msgMatrix,1,function(x) .encodePolynomial(x, genpoly, length, k))

  return(as.vector(retMatrix))

  #return(.encodePolynomial(msg,genpoly))

}

.BCHDecode = function(msg,params,visualize){

  if(length(params)>2)
    print("wrong params for bch, must be c(length,t)")

  length=params[1]
  t = params[2]
  m= ceiling(log2(length))


  if((length(msg)%%length)!= 0){
    msg = c(msg, rep(0,(length-(length(msg)%%length))));
  }


  k = .getGenPoly(m,length,t)$lm
 decMsg = .decodeBch(m,length,t,msg)

  msgMatrix = t(matrix(decMsg, nrow=length));
  retMatrix = apply(msgMatrix,1,function(x) x[(length-k+1):length])

  return(as.vector(retMatrix))

}

.encodePolynomial = function(data, g, length, k)
{

#   data_star = c(rep(0,length(g)-1),data)
#   res= polynomial(data_star)%%polynomial(g)
#   res = as.integer(res)
#  res = sapply(res ,function(x) x%%2)
#   res = as.integer(data_star - c(res,rep(0,length(data))))
#  res = sapply(res, function(x) x%%2)
#   return(res)

  dyn.load("src/bch.so")
  result=.C("encode_bch",length_r=as.integer(length), lm=as.integer(k), g=as.integer(g),bb=as.integer(rep(0,(length-k))), recd=as.integer(rep(0,length)), data=as.integer(data), PACKAGE = "bch")
  return(result$recd)

}

#'@useDynLib bch
.decodeBch = function(m,length,t,msg)
{
  dyn.load("src/bch.so")
  result=.C("decode_bch",m_r=as.integer(m),length_r=as.integer(length),t_r=as.integer(t),alpha_to=as.integer(rep(0,2**m)),index_of=as.integer(rep(0,2**m)),recd=as.integer(msg), nBlocks=as.integer(length(msg)/length), PACKAGE = "bch")
  return(result$recd)
}

#'@useDynLib bch
.getGenPoly= function(m, length, t)
{
  dyn.load("src/bch.so")
  result=.C("gen_poly",m_r=as.integer(m),length_r=as.integer(length),t_r=as.integer(t),g=as.integer(rep(0,2**m)),alpha_to=as.integer(rep(0,2**m)),index_of=as.integer(rep(0,2**m)),length_g=as.integer(0),lm=as.integer(0),PACKAGE = "bch")
  #return(c((result$g[1:(result$length_g+1)]),(result$lm)))
  return(result)
  #return(result$g[1:(result$length_g+1)])
}

