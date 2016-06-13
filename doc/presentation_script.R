block.coder <- BlockGenerateEncoderBCH(code.length = 15, code.t = 3)

# Faltungskodierer
#
#           ---->(+)-->(+)----->
#          |      ^     ^
#          |    __|_____|__
#   msg ---+-->|_____|_____|
#          |            |
#          |            v
#          ----------->(+)----->
#
conv.coder <- ConvGenerateEncoder(N = 2, M = 2, generators = c(7,5))

message <- sample(x = c(1,0), size = 100, replace = TRUE)

perm.vector = TurboGetPermutation(message.length = length(message),
                                  coder.info = conv.coder,
                                  type = "PRIMITIVE",
                                  args = list(root=2))

# encode
block.coded <- BlockEncode(message, block.encoder = block.coder)

conv.coded <- ConvEncode(message, conv.encoder = conv.coder)

turbo.coded <- TurboEncode(message, permutation.vector = perm.vector, coder.info = conv.coder)

# decode
block.decoded <- BlockDecode(block.coded, block.encoder = block.coder)

conv.decoded <- ConvDecodeSoft(conv.coded, conv.encoder = conv.coder)$output.hard

turbo.decoded <- TurboDecode(turbo.coded, permutation.vector = perm.vector,
                             iterations = 3, coder.info = conv.coder)$output.hard

# message = decoded ???
message
block.decoded
conv.decoded
turbo.decoded

sum(abs(message - block.decoded))
sum(abs(message - conv.decoded))
sum(abs(turbo.decoded - turbo.decoded))

# code with noise
block.noisy <- ApplyNoise(block.coded, SNR.db = 1, binary = TRUE)
block.messy <- ApplyNoise(block.coded, SNR.db = 0.1, binary = TRUE)

conv.noisy <- ApplyNoise(conv.coded, SNR.db = 1)
conv.messy <- ApplyNoise(conv.coded, SNR.db = 0.1)

turbo.noisy <- ApplyNoise(turbo.coded, SNR.db = 1)
turbo.messy <- ApplyNoise(turbo.coded, SNR.db = 0.1)

# show coded != noisy
block.coded
block.noisy
block.messy

conv.coded
conv.noisy
conv.messy

turbo.coded
turbo.noisy
turbo.messy

# decode noisy code
block.decoded <- BlockDecode(block.noisy, block.encoder = block.coder)

conv.decoded <- ConvDecodeSoft(conv.noisy, conv.encoder = conv.coder)$output.hard

turbo.decoded <- TurboDecode(turbo.noisy, permutation.vector = perm.vector,
                             iterations = 3, coder.info = conv.coder)$output.hard

# message = decoded ???
message
block.decoded
conv.decoded
turbo.decoded

sum(abs(message - block.decoded))
sum(abs(message - conv.decoded))
sum(abs(turbo.decoded - turbo.decoded))


##### visualizations #####

message = c(1, 0, 1, 0, 1)

block.coded <- BlockEncode(message, block.encoder = block.coder, visualize = TRUE)
block.decoded <- BlockDecode(block.coded, block.encoder = block.coder, visualize = TRUE)

# Evt. 2. Visualisierung, da man hier sieht was passiert wenn die Nachricht zu lang für einen Block ist
# und die Hamming Visualisierung sich von der BCH Visualisierung unterscheidet
#block.coder <- BlockGenerateEncoderHamming(code.length = 7, data.length = 4)
#block.coded <- BlockEncode(message, block.encoder = block.coder, visualize = TRUE)
#block.decoded <- BlockDecode(block.coded, block.encoder = block.coder, visualize = TRUE)

conv.coded <- ConvEncode(message, conv.encoder = conv.coder, visualize = TRUE)
#conv.decoded <- ConvDecodeSoft(conv.coded, conv.encoder = conv.coder, visualize = TRUE)

turbo.coded <- TurboEncode(message, permutation.vector = perm.vector,
                           coder.info = conv.coder, visualize = TRUE)
turbo.decoded <- TurboDecode(turbo.coded, permutation.vector = perm.vector,
                             iterations = 3, coder.info = conv.coder, visualize = TRUE)

##### simulation #####

df <- ChannelcodingSimulation(visualize = TRUE)

df
