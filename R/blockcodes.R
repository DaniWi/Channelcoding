.HammingEncode = function(msg, params, visualize)
{
  print("Hamming here!")
  return(c(msg,params))
}
.MatrixEncode = function(msg, params, visualize)
{
  print("Matrix here!")
  return(params%*%t(params))
}

