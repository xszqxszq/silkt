package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.entropy.RangeCoderState

internal class EncoderPayloadState {
    var frameCount: Int = 0
    var byteCount: Int = 0
    val rangeCoder: RangeCoderState = RangeCoderState()
}
