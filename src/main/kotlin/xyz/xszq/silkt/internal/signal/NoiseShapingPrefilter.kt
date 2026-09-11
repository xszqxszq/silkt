package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.MAX_SHAPE_LPC_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.model.SIG_TYPE_VOICED
import xyz.xszq.silkt.internal.quantization.HARMONIC_SHAPE_FIR_TAPS

// Ported from tsilk.
// Source: prefilter_FIX.ts

private const val LTP_SHAPING_BUFFER_MASK: Int = 511

@Suppress("SameParameterValue")
private fun warpedLpcAnalysisFilter(
    state: IntArray,
    residual: IntArray,
    residualOffset: Int,
    coefficientsQ13: IntArray,
    coefficientOffset: Int,
    input: IntArray,
    inputOffset: Int,
    warpingQ16: Int,
    length: Int,
    order: Int,
) {
    var firstTemporaryValue: Int
    var secondTemporaryValue: Int
    var sampleIndex = 0
    while (sampleIndex < length) {
        val inputSample = input[inputOffset + sampleIndex]
        secondTemporaryValue = smlawb(state[0], state[1], warpingQ16)
        state[0] = leftShift(inputSample, 14)
        firstTemporaryValue = smlawb(state[1], state[2] - secondTemporaryValue, warpingQ16)
        state[1] = secondTemporaryValue
        var accumulatorQ11 = smulwb(secondTemporaryValue, coefficientsQ13[coefficientOffset])

        var coefficientIndex = 2
        while (coefficientIndex < order) {
            secondTemporaryValue = smlawb(
                state[coefficientIndex],
                state[coefficientIndex + 1] - firstTemporaryValue,
                warpingQ16,
            )
            state[coefficientIndex] = firstTemporaryValue
            accumulatorQ11 = smlawb(
                accumulatorQ11,
                firstTemporaryValue,
                coefficientsQ13[coefficientOffset + coefficientIndex - 1],
            )
            firstTemporaryValue = smlawb(
                state[coefficientIndex + 1],
                state[coefficientIndex + 2] - secondTemporaryValue,
                warpingQ16,
            )
            state[coefficientIndex + 1] = secondTemporaryValue
            accumulatorQ11 = smlawb(
                accumulatorQ11,
                secondTemporaryValue,
                coefficientsQ13[coefficientOffset + coefficientIndex],
            )
            coefficientIndex += 2
        }

        state[order] = firstTemporaryValue
        accumulatorQ11 = smlawb(
            accumulatorQ11,
            firstTemporaryValue,
            coefficientsQ13[coefficientOffset + order - 1],
        )
        residual[residualOffset + sampleIndex] =
            sat16(inputSample - rshiftRound(accumulatorQ11, 11))
        sampleIndex++
    }
}

private fun prefilt(
    state: PrefilterState,
    shapedResidualQ12: IntArray,
    output: IntArray,
    outputOffset: Int,
    harmonicShapeFilterQ12: Int,
    tiltQ14: Int,
    lowFrequencyShapingQ14: Int,
    lag: Int,
    length: Int,
) {
    val ltpShapingBuffer = state.ltpShape
    var ltpShapingBufferIndex = state.ltpShapeWriteIndex
    var lowFrequencyARShapingQ12 = state.lowFrequencyShapingStateQ12
    var lowFrequencyMovingAverageShapingQ12 = state.lowFrequencyMovingAverageShapingQ12

    var sampleIndex = 0
    while (sampleIndex < length) {
        val ltpShapingQ12 = if (lag > 0) {
            val bufferIndex = lag + ltpShapingBufferIndex
            val filterOffset = HARMONIC_SHAPE_FIR_TAPS shr 1
            var shapedValue = smulbb(
                ltpShapingBuffer[bufferIndex - filterOffset - 1 and LTP_SHAPING_BUFFER_MASK],
                harmonicShapeFilterQ12,
            )
            shapedValue = smlabb(
                shapedValue,
                ltpShapingBuffer[bufferIndex - filterOffset and LTP_SHAPING_BUFFER_MASK],
                harmonicShapeFilterQ12 shr 16,
            )
            smlabb(
                shapedValue,
                ltpShapingBuffer[bufferIndex - filterOffset + 1 and LTP_SHAPING_BUFFER_MASK],
                harmonicShapeFilterQ12,
            )
        } else {
            0
        }

        val tiltQ10 = smulwb(lowFrequencyARShapingQ12, tiltQ14)
        val lowFrequencyQ10 = smlawb(
            smulwb(lowFrequencyARShapingQ12, lowFrequencyShapingQ14 shr 16),
            lowFrequencyMovingAverageShapingQ12,
            lowFrequencyShapingQ14,
        )
        lowFrequencyARShapingQ12 = sub32(
            shapedResidualQ12[sampleIndex],
            leftShift(tiltQ10, 2),
        )
        lowFrequencyMovingAverageShapingQ12 = sub32(
            lowFrequencyARShapingQ12,
            leftShift(lowFrequencyQ10, 2),
        )
        ltpShapingBufferIndex = ltpShapingBufferIndex - 1 and LTP_SHAPING_BUFFER_MASK
        ltpShapingBuffer[ltpShapingBufferIndex] =
            sat16(rshiftRound(lowFrequencyMovingAverageShapingQ12, 12))
        output[outputOffset + sampleIndex] =
            sat16(rshiftRound(sub32(lowFrequencyMovingAverageShapingQ12, ltpShapingQ12), 12))
        sampleIndex++
    }

    state.lowFrequencyShapingStateQ12 = lowFrequencyARShapingQ12
    state.lowFrequencyMovingAverageShapingQ12 = lowFrequencyMovingAverageShapingQ12
    state.ltpShapeWriteIndex = ltpShapingBufferIndex
}

internal fun prefilter(
    context: NoiseShapingContext,
    features: FrameFeatures,
    output: IntArray,
    input: IntArray,
    inputOffset: Int,
) {
    val state = context.prefilter
    val buffers = context.buffers
    val harmonicFilteredQ12 = buffers.harmonicFilteredQ12
    val shapedResidual = buffers.shapedResidual
    val subframeLength = context.subframeLength
    var inputPosition = inputOffset
    var outputOffset = 0
    var lag = state.previousLag

    for (subframeIndex in 0 until NB_SUBFR) {
        if (features.signalType == SIG_TYPE_VOICED) {
            lag = features.pitchL[subframeIndex]
        }
        val harmonicShapeGainQ12 = smulwb(
            features.harmShapeGainQ14[subframeIndex],
            16384 - features.harmBoostQ14[subframeIndex],
        )
        val harmonicShapeFilterQ12 = rightShift(harmonicShapeGainQ12, 2) or
                leftShift(rightShift(harmonicShapeGainQ12, 1), 16)
        val tiltQ14 = features.tiltQ14[subframeIndex]
        val lowFrequencyShapingQ14 = features.lowFrequencyShapingQ14[subframeIndex]

        warpedLpcAnalysisFilter(
            state.autoRegressiveShapingState,
            shapedResidual,
            0,
            features.ar1Q13,
            subframeIndex * MAX_SHAPE_LPC_ORDER,
            input,
            inputPosition,
            context.warpingQ16,
            subframeLength,
            context.shapingLpcOrder,
        )

        val gainQ12 = rshiftRound(features.gainsPreQ14[subframeIndex], 2)
        var temporaryValue = smlabb(
            fixConst(0.05, 26),
            features.harmBoostQ14[subframeIndex],
            harmonicShapeGainQ12,
        )
        temporaryValue = smlabb(
            temporaryValue,
            features.codingQualityQ14,
            fixConst(0.1, 12),
        )
        temporaryValue = smulwb(temporaryValue, -features.gainsPreQ14[subframeIndex])
        temporaryValue = rshiftRound(temporaryValue, 12)
        val harmonicCouplingQ12 = sat16(temporaryValue)

        harmonicFilteredQ12[0] = smlabb(
            smulbb(shapedResidual[0], gainQ12),
            state.harmonicHighPassState,
            harmonicCouplingQ12,
        )
        for (sampleIndex in 1 until subframeLength) {
            harmonicFilteredQ12[sampleIndex] = smlabb(
                smulbb(shapedResidual[sampleIndex], gainQ12),
                shapedResidual[sampleIndex - 1],
                harmonicCouplingQ12,
            )
        }
        state.harmonicHighPassState = shapedResidual[subframeLength - 1]

        prefilt(
            state,
            harmonicFilteredQ12,
            output,
            outputOffset,
            harmonicShapeFilterQ12,
            tiltQ14,
            lowFrequencyShapingQ14,
            lag,
            subframeLength,
        )
        inputPosition += subframeLength
        outputOffset += subframeLength
    }
    state.previousLag = features.pitchL[NB_SUBFR - 1]
}
