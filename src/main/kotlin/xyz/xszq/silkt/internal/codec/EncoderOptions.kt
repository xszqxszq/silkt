package xyz.xszq.silkt.internal.codec

internal class EncoderOptions(
    var apiSampleRate: Int,
    var maxInternalSampleRate: Int? = null,
    var packetSize: Int? = null,
    var bitRate: Int? = null,
    var packetLossPercentage: Int? = null,
    var complexity: Int? = null,
    var useInBandFec: Boolean? = null,
    var useDtx: Boolean? = null
)
