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

message <- rep(c(1,0,1,1,1,0,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,0,1,0,0), times = 4)
message

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
message - block.decoded
message - conv.decoded
message - turbo.decoded

# code with noise
block.noisy <- ApplyNoise(block.coded, SNR.db = 1, binary = TRUE)
block.messy <- ApplyNoise(block.coded, SNR.db = 0.1, binary = TRUE)

conv.noisy <- ApplyNoise(conv.coded, SNR.db = 1)
conv.messy <- ApplyNoise(conv.coded, SNR.db = 0.1)

turbo.noisy <- ApplyNoise(turbo.coded, SNR.db = 1)
turbo.messy <- ApplyNoise(turbo.coded, SNR.db = 0.1)

# show coded != noisy
block.coded - block.noisy
block.coded - block.messy

conv.coded
round(conv.noisy, digits = 2)
round(conv.messy, digits = 2)

turbo.coded
round(turbo.noisy, digits = 2)
round(turbo.messy, digits = 2)

# decode noisy code
# noisy & messy!!!
block.noisy.decoded <- BlockDecode(block.noisy, block.encoder = block.coder)
block.messy.decoded <- BlockDecode(block.messy, block.encoder = block.coder)

conv.noisy.decoded <- ConvDecodeSoft(conv.noisy, conv.encoder = conv.coder)$output.hard
conv.messy.decoded <- ConvDecodeSoft(conv.messy, conv.encoder = conv.coder)$output.hard

turbo.noisy.decoded <- TurboDecode(turbo.noisy, permutation.vector = perm.vector,
                             iterations = 3, coder.info = conv.coder)$output.hard
turbo.noisy.decoded <- TurboDecode(turbo.messy, permutation.vector = perm.vector,
                             iterations = 3, coder.info = conv.coder)$output.hard

# message = decoded ???
message - block.noisy.decoded
message - block.messy.decoded

message - conv.noisy.decoded
message - conv.messy.decoded

message - turbo.noisy.decoded
message - turbo.messy.decoded

##### visualizations #####

message = c(1, 0, 1, 0, 1)

block.coded <- BlockEncode(message, block.encoder = block.coder, visualize = TRUE)
block.decoded <- BlockDecode(block.coded, block.encoder = block.coder, visualize = TRUE)

conv.coded <- ConvEncode(message, conv.encoder = conv.coder, visualize = TRUE)
# conv.decoded <- ConvDecodeSoft(conv.coded, conv.encoder = conv.coder)

turbo.coded <- TurboEncode(message, permutation.vector = perm.vector,
                           coder.info = conv.coder)
turbo.decoded <- TurboDecode(turbo.coded, permutation.vector = perm.vector,
                             iterations = 3, coder.info = conv.coder, visualize = TRUE)

##### simulation #####

df <- ChannelcodingSimulation(visualize = TRUE)
df
