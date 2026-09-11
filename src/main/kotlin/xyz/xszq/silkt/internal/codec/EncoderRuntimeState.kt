package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.util.RefInt

internal class EncoderRuntimeState {
    var frameCounter: Int = 0
    var voiceActivityFlag: Boolean = false
    var inDtx: Boolean = false
    var silentFrameCount: Int = 0
    val previousTypeOffset: RefInt = RefInt()
    var firstFrameAfterReset: Boolean = false
    var speechActivityQ8: Int = 0
    var channelBufferedMs: Int = 0
}
