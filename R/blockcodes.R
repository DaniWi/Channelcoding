.HammingEncode = function(msg, params, visualize)
{
  print("Hamming here!")



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

