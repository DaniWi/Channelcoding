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

systematic.coder <- ConvGenerateRscEncoder(2,2,c(5,7))

message <- rep(c(1,0,1,1,1,0,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,0,1,0,0), times = 4)
message

perm.vector = TurboGetPermutation(message.length = length(message),
                                  coder.info = conv.coder,
                                  type = "PRIMITIVE",
                                  list(root=2),
                                  visualize = TRUE)

# encode
block.coded <- BlockEncode(message, block.encoder = block.coder)
block.coded

conv.coded <- ConvEncode(message, conv.encoder = conv.coder)
conv.coded

turbo.coded <- TurboEncode(message, permutation.vector = perm.vector, coder.info = systematic.coder)

# decode
block.decoded <- BlockDecode(block.coded, block.encoder = block.coder)

conv.decoded <- ConvDecodeSoft(conv.coded, conv.encoder = conv.coder)$output.hard

turbo.decoded <- TurboDecode(turbo.coded, permutation.vector = perm.vector,
                             iterations = 3, coder.info = systematic.coder)$output.hard

# message = decoded ???
message - block.decoded
message - conv.decoded
message - turbo.decoded

# code with noise
Block.noisy <- ApplyNoise(block.coded, SNR.db = 2, binary = TRUE)
Block.messy <- ApplyNoise(block.coded, SNR.db = 0.1, binary = TRUE)

Conv.noisy <- ApplyNoise(conv.coded, SNR.db = 2)
Conv.messy <- ApplyNoise(conv.coded, SNR.db = 0.1)

Turbo.noisy <- ApplyNoise(turbo.coded, SNR.db = 2)
Turbo.messy <- ApplyNoise(turbo.coded, SNR.db = 0.1)

# show coded != noisy
block.coded - block.noisy
block.coded - block.messy

round(conv.noisy, digits = 2)
#round(conv.messy, digits = 2)

round(turbo.noisy, digits = 2)
#round(turbo.messy, digits = 2)

# decode noisy code
block.noisy.decoded <- BlockDecode(block.noisy, block.encoder = block.coder)
block.messy.decoded <- BlockDecode(block.messy, block.encoder = block.coder)

conv.noisy.decoded <- ConvDecodeSoft(conv.noisy, conv.encoder = conv.coder)$output.hard
conv.messy.decoded <- ConvDecodeSoft(conv.messy, conv.encoder = conv.coder)$output.hard

turbo.noisy.decoded <- TurboDecode(turbo.noisy, permutation.vector = perm.vector,
                             iterations = 10, coder.info = systematic.coder)$output.hard
turbo.messy.decoded <- TurboDecode(turbo.messy, permutation.vector = perm.vector,
                             iterations = 10, coder.info = systematic.coder)$output.hard

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
#block.decoded <- BlockDecode(block.coded, block.encoder = block.coder)

conv.coded <- ConvEncode(message, conv.encoder = conv.coder, visualize = TRUE)

perm.vector = TurboGetPermutation(message.length = length(message),
                                  coder.info = conv.coder,
                                  type = "RANDOM")

turbo.coded <- TurboEncode(message, permutation.vector = perm.vector,
                           coder.info = systematic.coder)
turbo.noisy.vis <- ApplyNoise(turbo.coded, SNR.db = 3)
turbo.decoded <- TurboDecode(turbo.noisy.vis, permutation.vector = perm.vector,
                             iterations = 3, coder.info = systematic.coder, visualize = TRUE)

##### simulation #####

df <- ChannelcodingSimulation(visualize = TRUE)
df
