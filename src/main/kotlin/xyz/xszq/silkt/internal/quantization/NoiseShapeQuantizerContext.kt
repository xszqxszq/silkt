package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.util.RefInt

internal class NoiseShapeQuantizerContext {
    var frameLength: Int = 0
    var subframeLength: Int = 0
    var predictionLpcOrder: Int = 16
    var shapingLpcOrder: Int = 16
    var delayedDecisionCount: Int = 1
    var warpingQ16: Int = 0
    val primaryState: NoiseShapeQuantizerState = NoiseShapeQuantizerState()
    val lbrrState: NoiseShapeQuantizerState = NoiseShapeQuantizerState()
    val scaledLtpExcitationQ16 = IntArray(2 * MAX_FRAME_LENGTH)
    val ltpExcitation = IntArray(2 * MAX_FRAME_LENGTH)
    val scaledInputQ10 = IntArray(MAX_FRAME_LENGTH / NB_SUBFR)
    val lpcFilterState = IntArray(MAX_LPC_ORDER)
    val delayedDecisions = Array(MAX_DELAYED_DECISION_STATES) { createDelayedDecisionState() }
    val sampleStates = Array(MAX_DELAYED_DECISION_STATES) { Array(2) { SampleState() } }
    val sampleHistoryIndex = RefInt()

    fun configure(
        frameLength: Int,
        subframeLength: Int,
        predictionLpcOrder: Int,
        shapingLpcOrder: Int,
        delayedDecisionCount: Int,
        warpingQ16: Int,
    ) {
        this.frameLength = frameLength
        this.subframeLength = subframeLength
        this.predictionLpcOrder = predictionLpcOrder
        this.shapingLpcOrder = shapingLpcOrder
        this.delayedDecisionCount = delayedDecisionCount
        this.warpingQ16 = warpingQ16
    }
}
