package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.model.SIG_TYPE_VOICED
import xyz.xszq.silkt.internal.signal.FrameFeatures
import xyz.xszq.silkt.internal.signal.ShapeState
import xyz.xszq.silkt.internal.signal.sigmoidQ15
import xyz.xszq.silkt.internal.tables.ExcitationTables.quantizationOffsetsQ10
import xyz.xszq.silkt.internal.util.RefInt

private const val MIN_QUANTIZED_GAIN_DB = 6
private const val MAX_QUANTIZED_GAIN_DB = 86
private const val MAX_DELTA_GAIN_QUANT = 40
private const val MIN_DELTA_GAIN_QUANT = -4
private const val TUNING_LAMBDA_OFFSET = 1.2
private const val TUNING_LAMBDA_SPEECH_ACT = -0.3
private const val TUNING_LAMBDA_DELAYED_DECISIONS = -0.05
private const val TUNING_LAMBDA_INPUT_QUALITY = -0.2
private const val TUNING_LAMBDA_CODING_QUALITY = -0.1
private const val TUNING_LAMBDA_QUANT_OFFSET = 1.5
private const val GAIN_QUANTIZATION_OFFSET_Q7 = 2176
private const val QUANTIZED_GAIN_RANGE_Q7 = ((MAX_QUANTIZED_GAIN_DB - MIN_QUANTIZED_GAIN_DB) * 128) / 6

internal const val GAIN_LEVEL_COUNT = 64

private const val GAIN_SCALE_Q16 =
    (65536 * (GAIN_LEVEL_COUNT - 1)) / QUANTIZED_GAIN_RANGE_Q7

private const val INVERSE_GAIN_SCALE_Q16 =
    (65536 * QUANTIZED_GAIN_RANGE_Q7) / (GAIN_LEVEL_COUNT - 1)

// Ported from tsilk.
// Source: process_gains_FIX.ts

internal fun processGains(
    shapeState: ShapeState,
    subframeLength: Int,
    payloadFrameCount: Int,
    delayedDecisionCount: Int,
    speechActivityQ8: Int,
    features: FrameFeatures,
) {
    var subframeIndex = 0
    var smoothingGainQ16: Int
    val inverseMaximumSquareValueQ16 = div3216(
        log2Lin(
            smulwb(
                fixConst(70.0, 7) - features.currentSnrDbQ7,
                fixConst(0.33, 16),
            ),
        ),
        subframeLength,
    )

    if (features.signalType == SIG_TYPE_VOICED) {
        smoothingGainQ16 = -sigmoidQ15(rshiftRound(features.ltPredCodGainQ7 - fixConst(12.0, 7), 4))
        while (subframeIndex < NB_SUBFR) {
            features.gainsQ16[subframeIndex] = smlawb(
                features.gainsQ16[subframeIndex],
                features.gainsQ16[subframeIndex],
                smoothingGainQ16,
            )
            subframeIndex++
        }
    }

    subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val residualEnergy = features.residualEnergy[subframeIndex]
        var residualEnergyPart = smulww(residualEnergy, inverseMaximumSquareValueQ16)
        val residualEnergyShift = features.residualEnergyShift[subframeIndex]

        if (residualEnergyShift > 0) {
            residualEnergyPart =
                if (residualEnergyShift < 32) rshiftRound(residualEnergyPart, residualEnergyShift) else 0
        } else if (residualEnergyShift != 0) {
            residualEnergyPart = if (residualEnergyPart > rightShift(Int.MAX_VALUE, -residualEnergyShift)) {
                Int.MAX_VALUE
            } else {
                leftShift(residualEnergyPart, -residualEnergyShift)
            }
        }

        var gain = features.gainsQ16[subframeIndex]
        val gainSquared = residualEnergyPart + smmul(gain, gain)
        if (gainSquared < 32767) {
            gain = sqrtApprox(smlaww(leftShift(residualEnergyPart, 16), gain, gain))
            features.gainsQ16[subframeIndex] = lshiftSat32(gain, 8)
        } else {
            gain = sqrtApprox(gainSquared)
            features.gainsQ16[subframeIndex] = lshiftSat32(gain, 16)
        }
        subframeIndex++
    }

    val previousGainIndex = RefInt(shapeState.lastGainIndex)
    quantizeGains(
        features.gainsIndices,
        features.gainsQ16,
        previousGainIndex,
        if (payloadFrameCount > 0) 1 else 0,
    )
    shapeState.lastGainIndex = previousGainIndex.value

    if (features.signalType == SIG_TYPE_VOICED) {
        features.quantOffsetType =
            if (features.ltPredCodGainQ7 + rightShift(features.inputTiltQ15, 8) > fixConst(1.0, 7)) 0 else 1
    }

    val quantizationOffsetQ10 = quantizationOffsetsQ10[features.signalType][features.quantOffsetType]
    features.lambdaQ10 =
        fixConst(TUNING_LAMBDA_OFFSET, 10) +
                smulbb(fixConst(TUNING_LAMBDA_DELAYED_DECISIONS, 10), delayedDecisionCount) +
                smulwb(fixConst(TUNING_LAMBDA_SPEECH_ACT, 18), speechActivityQ8) +
                smulwb(fixConst(TUNING_LAMBDA_INPUT_QUALITY, 12), features.inputQualityQ14) +
                smulwb(fixConst(TUNING_LAMBDA_CODING_QUALITY, 12), features.codingQualityQ14) +
                smulwb(fixConst(TUNING_LAMBDA_QUANT_OFFSET, 16), quantizationOffsetQ10)
}

// Ported from tsilk.
// Source: gain_quant.ts

private fun quantizeGains(
    gainIndices: IntArray,
    gainsQ16: IntArray,
    previousGainIndex: RefInt,
    conditionalCoding: Int,
) {
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        gainIndices[subframeIndex] = smulwb(
            GAIN_SCALE_Q16,
            lin2Log(gainsQ16[subframeIndex]) - GAIN_QUANTIZATION_OFFSET_Q7,
        )
        if (gainIndices[subframeIndex] < previousGainIndex.value) {
            gainIndices[subframeIndex]++
        }

        if (subframeIndex == 0 && conditionalCoding == 0) {
            gainIndices[subframeIndex] = limitInt(gainIndices[subframeIndex], 0, GAIN_LEVEL_COUNT - 1)
            gainIndices[subframeIndex] =
                maxInt(gainIndices[subframeIndex], previousGainIndex.value + MIN_DELTA_GAIN_QUANT)
            previousGainIndex.value = gainIndices[subframeIndex]
        } else {
            gainIndices[subframeIndex] = limitInt(
                gainIndices[subframeIndex] - previousGainIndex.value,
                MIN_DELTA_GAIN_QUANT,
                MAX_DELTA_GAIN_QUANT,
            )
            previousGainIndex.value += gainIndices[subframeIndex]
            gainIndices[subframeIndex] -= MIN_DELTA_GAIN_QUANT
        }

        gainsQ16[subframeIndex] = log2Lin(
            min32(
                toInt32(smulwb(INVERSE_GAIN_SCALE_Q16, previousGainIndex.value) + GAIN_QUANTIZATION_OFFSET_Q7),
                3967,
            ),
        )
        subframeIndex++
    }
}

internal fun dequantizeGains(
    gainsQ16: IntArray,
    gainIndices: IntArray,
    previousGainIndex: RefInt,
    conditionalCoding: Int,
) {
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        previousGainIndex.value =
            if (subframeIndex == 0 && conditionalCoding == 0) {
                gainIndices[subframeIndex]
            } else {
                previousGainIndex.value + gainIndices[subframeIndex] + MIN_DELTA_GAIN_QUANT
            }

        gainsQ16[subframeIndex] = log2Lin(
            min32(
                toInt32(smulwb(INVERSE_GAIN_SCALE_Q16, previousGainIndex.value) + GAIN_QUANTIZATION_OFFSET_Q7),
                3967,
            ),
        )
        subframeIndex++
    }
}
