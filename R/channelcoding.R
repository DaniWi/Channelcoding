encode = function(msg, type, params, visualize)
{
  switch(type,
         HAMMING={
           return(.HammingEncode(msg,params,visualize))
         },
         MATRIX={
           return(.MatrixEncode(msg,params,visualize))
         },
         {
           print("No or Wrong Type selected")
           return()
         })

}

decode = function(msg, type, params, visualize)
{
  return(msg)
}

applyNoise = function(msg, params, visualize)
{
  return(msg)
}

