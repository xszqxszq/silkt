package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.model.NB_SUBFR

internal class EncoderFrameGeometry {
    var samplingRateKHz: Int = 16
    var frameLength: Int = 0
    var subframeLength: Int = 0
    var pitchLookahead: Int = 0
    var shapeLookahead: Int = 0

    fun configure(samplingRateKHz: Int, pitchLookaheadMs: Int, shapeLookaheadMs: Int) {
        this.samplingRateKHz = samplingRateKHz
        frameLength = FRAME_LENGTH_MS * samplingRateKHz
        subframeLength = frameLength / NB_SUBFR
        pitchLookahead = pitchLookaheadMs * samplingRateKHz
        shapeLookahead = shapeLookaheadMs * samplingRateKHz
    }
}
