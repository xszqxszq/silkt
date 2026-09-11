package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.quantization.NlsfCodebook

internal class EncoderConfig {
    var complexity: Int = 2
    var targetRateBps: Int = 25000
    var packetSizeMs: Int = FRAME_LENGTH_MS
    var packetLossPercentage: Int = 0
    var useInBandFec: Boolean = false
    var useDtx: Boolean = false
    var predictionLpcOrder: Int = 16
    var nlsfCodebooks: Array<NlsfCodebook> = arrayOf()
    var delayedDecisionCount: Int = 1
    var lbrrEnabled: Boolean = false
    var lbrrGainStepCount: Int = 0
    var inBandFecSnrCompensationQ8: Int = 0
    var targetSnrDbQ7: Int = 0
}
