package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.LTP_ORDER
import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.model.MAX_SHAPE_LPC_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR

internal const val DECISION_DELAY = 32
internal const val LPC_HISTORY_BUFFER_LENGTH = DECISION_DELAY
internal const val HARMONIC_SHAPE_FIR_TAPS = 3
internal const val MAX_DELAYED_DECISION_STATES = 4
internal const val DECISION_DELAY_MASK = DECISION_DELAY - 1

// Ported from tsilk.
// Source: Nsq_del_dec.ts

internal class NoiseShapeQuantizerState(
    var quantizedSignalQ10: IntArray = IntArray(2 * MAX_FRAME_LENGTH),
    var ltpShapeQ10: IntArray = IntArray(2 * MAX_FRAME_LENGTH),
    var lpcStateQ14: IntArray = IntArray(LPC_HISTORY_BUFFER_LENGTH + MAX_FRAME_LENGTH / NB_SUBFR),
    var autoregressiveStateQ14: IntArray = IntArray(MAX_SHAPE_LPC_ORDER),
    var lowFrequencyShapingStateQ12: Int = 0,
    var previousLag: Int = 100,
    var ltpExcitationWriteIndex: Int = 0,
    var ltpShapeWriteIndex: Int = 0,
    var randomSeed: Int = 0,
    var previousInverseGainQ16: Int = 0,
    var shouldRewhiten: Boolean = false,
)

internal fun createDelayedDecisionState(): DelayedDecisionState =
    DelayedDecisionState(
        seedHistory = IntArray(DECISION_DELAY),
        pulseQ10 = IntArray(DECISION_DELAY),
        quantizedSignalQ10 = IntArray(DECISION_DELAY),
        lpcExcitationQ16 = IntArray(DECISION_DELAY),
        ltpShapeQ10 = IntArray(DECISION_DELAY),
        gainQ16 = IntArray(DECISION_DELAY),
        autoregressiveStateQ14 = IntArray(MAX_SHAPE_LPC_ORDER),
        lpcStateQ14 = IntArray(LPC_HISTORY_BUFFER_LENGTH + MAX_FRAME_LENGTH / NB_SUBFR),
    )

internal fun selectNoiseShapeQuantizerState(
    context: NoiseShapeQuantizerContext,
    useLbrr: Boolean = false,
): NoiseShapeQuantizerState {
    val state = if (useLbrr) context.lbrrState else context.primaryState
    if (state.previousInverseGainQ16 == 0) {
        state.previousInverseGainQ16 = 65536
    }
    return state
}

internal fun copyFinalNoiseShapeQuantizerState(state: NoiseShapeQuantizerState, frameLength: Int) {
    state.quantizedSignalQ10.copyInto(state.quantizedSignalQ10, 0, frameLength, 2 * frameLength)
    state.ltpShapeQ10.copyInto(state.ltpShapeQ10, 0, frameLength, 2 * frameLength)
}

private fun rewhitenLtpExcitation(
    state: NoiseShapeQuantizerState,
    ltpExcitation: IntArray,
    scaledLtpExcitationQ16: IntArray,
    lag: Int,
    inverseGainQ16: Int,
    subframeIndex: Int,
    ltpScaleQ14: Int,
) {
    var inverseGainQ32 = leftShift(inverseGainQ16, 16)
    if (subframeIndex == 0) {
        inverseGainQ32 = leftShift(smulwb(inverseGainQ32, ltpScaleQ14), 2)
    }
    var sampleIndex = state.ltpExcitationWriteIndex - lag - (LTP_ORDER shr 1)
    while (sampleIndex < state.ltpExcitationWriteIndex) {
        scaledLtpExcitationQ16[sampleIndex] = smulwb(inverseGainQ32, ltpExcitation[sampleIndex])
        sampleIndex++
    }
}

private fun rescaleLtpShapeHistory(
    state: NoiseShapeQuantizerState,
    subframeLength: Int,
    gainAdjustmentQ16: Int,
) {
    var shapeIndex = state.ltpShapeWriteIndex - subframeLength * NB_SUBFR
    while (shapeIndex < state.ltpShapeWriteIndex) {
        state.ltpShapeQ10[shapeIndex] = smulww(gainAdjustmentQ16, state.ltpShapeQ10[shapeIndex])
        shapeIndex++
    }
}

private fun rescaleLtpExcitationHistory(
    state: NoiseShapeQuantizerState,
    scaledLtpExcitationQ16: IntArray,
    lag: Int,
    gainAdjustmentQ16: Int,
) {
    var sampleIndex = state.ltpExcitationWriteIndex - lag - (LTP_ORDER shr 1)
    while (sampleIndex < state.ltpExcitationWriteIndex) {
        scaledLtpExcitationQ16[sampleIndex] =
            smulww(gainAdjustmentQ16, scaledLtpExcitationQ16[sampleIndex])
        sampleIndex++
    }
}

@Suppress("DuplicatedCode")
private fun rescaleQuantizerStateHistory(
    state: NoiseShapeQuantizerState,
    gainAdjustmentQ16: Int,
) {
    state.lowFrequencyShapingStateQ12 = smulww(gainAdjustmentQ16, state.lowFrequencyShapingStateQ12)
    var historyIndex = 0
    while (historyIndex < LPC_HISTORY_BUFFER_LENGTH) {
        state.lpcStateQ14[historyIndex] = smulww(gainAdjustmentQ16, state.lpcStateQ14[historyIndex])
        historyIndex++
    }

    historyIndex = 0
    while (historyIndex < MAX_SHAPE_LPC_ORDER) {
        state.autoregressiveStateQ14[historyIndex] =
            smulww(gainAdjustmentQ16, state.autoregressiveStateQ14[historyIndex])
        historyIndex++
    }
}

@Suppress("DuplicatedCode")
private fun rescaleDelayedDecisionHistory(
    delayedDecision: DelayedDecisionState,
    gainAdjustmentQ16: Int,
) {
    delayedDecision.lowFrequencyShapingStateQ12 =
        smulww(gainAdjustmentQ16, delayedDecision.lowFrequencyShapingStateQ12)

    var historyIndex = 0
    while (historyIndex < LPC_HISTORY_BUFFER_LENGTH) {
        delayedDecision.lpcStateQ14[historyIndex] =
            smulww(gainAdjustmentQ16, delayedDecision.lpcStateQ14[historyIndex])
        historyIndex++
    }

    historyIndex = 0
    while (historyIndex < MAX_SHAPE_LPC_ORDER) {
        delayedDecision.autoregressiveStateQ14[historyIndex] =
            smulww(gainAdjustmentQ16, delayedDecision.autoregressiveStateQ14[historyIndex])
        historyIndex++
    }

    historyIndex = 0
    while (historyIndex < DECISION_DELAY) {
        delayedDecision.lpcExcitationQ16[historyIndex] =
            smulww(gainAdjustmentQ16, delayedDecision.lpcExcitationQ16[historyIndex])
        delayedDecision.ltpShapeQ10[historyIndex] =
            smulww(gainAdjustmentQ16, delayedDecision.ltpShapeQ10[historyIndex])
        historyIndex++
    }
}

private fun scaleQuantizerInput(
    input: IntArray,
    inputOffset: Int,
    scaledInputQ10: IntArray,
    subframeLength: Int,
    inverseGainQ16: Int,
) {
    var sampleIndex = 0
    while (sampleIndex < subframeLength) {
        scaledInputQ10[sampleIndex] =
            rightShift(smulbb(input[inputOffset + sampleIndex], inverseGainQ16), 6)
        sampleIndex++
    }
}

internal fun copyDelayedDecisionState(
    destination: DelayedDecisionState,
    source: DelayedDecisionState,
    lpcHistoryOffset: Int,
) {
    source.seedHistory.copyInto(destination.seedHistory)
    source.pulseQ10.copyInto(destination.pulseQ10)
    source.lpcExcitationQ16.copyInto(destination.lpcExcitationQ16)
    source.ltpShapeQ10.copyInto(destination.ltpShapeQ10)
    source.quantizedSignalQ10.copyInto(destination.quantizedSignalQ10)
    source.autoregressiveStateQ14.copyInto(destination.autoregressiveStateQ14)
    source.lpcStateQ14.copyInto(
        destination.lpcStateQ14,
        lpcHistoryOffset,
        lpcHistoryOffset,
        lpcHistoryOffset + LPC_HISTORY_BUFFER_LENGTH,
    )
    destination.lowFrequencyShapingStateQ12 = source.lowFrequencyShapingStateQ12
    destination.seed = source.seed
    destination.frameStartSeed = source.frameStartSeed
    destination.rateDistortionQ10 = source.rateDistortionQ10
}

internal fun scaleNoiseShapeQuantizerState(
    state: NoiseShapeQuantizerState,
    input: IntArray,
    inputOffset: Int,
    scaledInputQ10: IntArray,
    subframeLength: Int,
    ltpExcitation: IntArray,
    scaledLtpExcitationQ16: IntArray,
    subframeIndex: Int,
    ltpScaleQ14: Int,
    gainsQ16: IntArray,
    pitchLags: IntArray,
): Int {
    val inverseGainQ16 = min(inverse32VarQ(max(gainsQ16[subframeIndex], 1), 32), 32767)
    val lag = pitchLags[subframeIndex]

    if (state.shouldRewhiten) {
        rewhitenLtpExcitation(
            state,
            ltpExcitation,
            scaledLtpExcitationQ16,
            lag,
            inverseGainQ16,
            subframeIndex,
            ltpScaleQ14
        )
    }

    if (inverseGainQ16 != state.previousInverseGainQ16) {
        val gainAdjustmentQ16 = div32VarQ(inverseGainQ16, state.previousInverseGainQ16, 16)
        rescaleLtpShapeHistory(state, subframeLength, gainAdjustmentQ16)

        if (!state.shouldRewhiten) {
            rescaleLtpExcitationHistory(state, scaledLtpExcitationQ16, lag, gainAdjustmentQ16)
        }

        rescaleQuantizerStateHistory(state, gainAdjustmentQ16)
    }

    scaleQuantizerInput(input, inputOffset, scaledInputQ10, subframeLength, inverseGainQ16)
    state.previousInverseGainQ16 = inverseGainQ16
    return inverseGainQ16
}

internal fun scaleDelayedDecisionStates(
    state: NoiseShapeQuantizerState,
    delayedDecisions: Array<DelayedDecisionState>,
    input: IntArray,
    inputOffset: Int,
    scaledInputQ10: IntArray,
    subframeLength: Int,
    ltpExcitation: IntArray,
    scaledLtpExcitationQ16: IntArray,
    subframeIndex: Int,
    delayedDecisionCount: Int,
    ltpScaleQ14: Int,
    gainsQ16: IntArray,
    pitchLags: IntArray,
): Int {
    val inverseGainQ16 = min(inverse32VarQ(max(gainsQ16[subframeIndex], 1), 32), 32767)
    val lag = pitchLags[subframeIndex]

    if (state.shouldRewhiten) {
        rewhitenLtpExcitation(
            state,
            ltpExcitation,
            scaledLtpExcitationQ16,
            lag,
            inverseGainQ16,
            subframeIndex,
            ltpScaleQ14
        )
    }

    if (inverseGainQ16 != state.previousInverseGainQ16) {
        val gainAdjustmentQ16 = div32VarQ(inverseGainQ16, state.previousInverseGainQ16, 16)
        rescaleLtpShapeHistory(state, subframeLength, gainAdjustmentQ16)

        if (!state.shouldRewhiten) {
            rescaleLtpExcitationHistory(state, scaledLtpExcitationQ16, lag, gainAdjustmentQ16)
        }

        var decisionIndex = 0
        while (decisionIndex < delayedDecisionCount) {
            val delayedDecision = delayedDecisions[decisionIndex]
            rescaleDelayedDecisionHistory(delayedDecision, gainAdjustmentQ16)
            decisionIndex++
        }
    }

    scaleQuantizerInput(input, inputOffset, scaledInputQ10, subframeLength, inverseGainQ16)
    state.previousInverseGainQ16 = inverseGainQ16
    return inverseGainQ16
}
