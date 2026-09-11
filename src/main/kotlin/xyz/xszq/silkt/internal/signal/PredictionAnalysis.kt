package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.model.SIG_TYPE_VOICED
import xyz.xszq.silkt.internal.quantization.computeLaroiaNlsfWeights
import xyz.xszq.silkt.internal.quantization.convertNlsfToStablePrediction
import xyz.xszq.silkt.internal.quantization.encodeNlsfMsvq
import xyz.xszq.silkt.internal.quantization.interpolateValues
import xyz.xszq.silkt.internal.util.RefInt

internal class PredictionState(
    var prevNlsfQ15: IntArray = IntArray(16),
)

// Ported from tsilk.
// Source: find_pred_coefs_FIX.ts

private fun ensurePreviousNlsf(context: PredictionAnalysisContext, lpcOrder: Int): IntArray {
    if (context.state.prevNlsfQ15.size != lpcOrder) {
        context.state.prevNlsfQ15 = IntArray(lpcOrder)
    }
    return context.state.prevNlsfQ15
}

internal fun findPredictionCoefficients(
    context: PredictionAnalysisContext,
    features: FrameFeatures,
    residualPitch: IntArray,
    inputBuffer: IntArray,
    firstFrameAfterReset: Boolean,
    speechActivityQ8: Int,
    packetLossPercentage: Int,
    framesInPayloadBuffer: Int,
    packetSizeMilliseconds: Int,
) {
    val workspace = context.reusableWorkspace
    val lpcOrder = if (context.predictionLpcOrder != 0) context.predictionLpcOrder else MAX_LPC_ORDER
    val subframeLength = context.subframeLength
    val inverseGainsQ16 = workspace.inverseGainsQ16
    val localGains = workspace.localGains

    features.predCoefQ12.fill(0)
    features.ltpCoefQ14.fill(0)
    features.ltpScaleQ14 = 0
    features.ltpIndex.fill(0)
    features.nlsfIndices.fill(0)

    prepareSubframeGains(features, inverseGainsQ16, localGains, workspace.subframeWeightsQ15)
    val lpcInput = workspace.lpcInput
    buildPredictionInput(
        context,
        features,
        workspace,
        residualPitch,
        inputBuffer,
        inverseGainsQ16,
        lpcOrder,
        subframeLength,
        packetLossPercentage,
        framesInPayloadBuffer,
        packetSizeMilliseconds,
    )

    val previousNlsfQ15 = ensurePreviousNlsf(context, lpcOrder)
    val nlsfQ15 = workspace.nlsfQ15
    val interpolationIndex = analyzeNlsf(
        context,
        features,
        workspace,
        lpcInput,
        previousNlsfQ15,
        nlsfQ15,
        firstFrameAfterReset,
        speechActivityQ8,
        lpcOrder,
        subframeLength,
    )
    buildPredictionCoefficients(
        features,
        workspace,
        previousNlsfQ15,
        nlsfQ15,
        interpolationIndex,
        lpcOrder
    )

    residualEnergy(
        workspace,
        features.residualEnergy,
        features.residualEnergyShift,
        lpcInput,
        0,
        localGains,
        subframeLength,
        lpcOrder
    )
    nlsfQ15.copyInto(previousNlsfQ15)
    restoreZeroGains(features, subframeLength)
}

private fun prepareSubframeGains(
    features: FrameFeatures,
    inverseGainsQ16: IntArray,
    localGains: IntArray,
    subframeWeightsQ15: IntArray
) {
    var minimumGainQ16 = Int.MAX_VALUE shr 6
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        minimumGainQ16 = min(minimumGainQ16, max(1, features.gainsQ16[subframeIndex]))
        subframeIndex++
    }

    subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val gainQ16 = max(1, features.gainsQ16[subframeIndex])
        var inverseGainQ16 = div32VarQ(minimumGainQ16, gainQ16, 14)
        inverseGainQ16 = max(inverseGainQ16, 363)
        inverseGainsQ16[subframeIndex] = inverseGainQ16
        subframeWeightsQ15[subframeIndex] = rightShift(
            smulwb(inverseGainQ16, inverseGainQ16),
            1
        )
        localGains[subframeIndex] = div32(1 shl 16, inverseGainQ16)
        subframeIndex++
    }
}

private fun buildPredictionInput(
    context: PredictionAnalysisContext,
    features: FrameFeatures,
    workspace: PredictionAnalysisWorkspace,
    residualPitch: IntArray,
    inputBuffer: IntArray,
    inverseGainsQ16: IntArray,
    lpcOrder: Int,
    subframeLength: Int,
    packetLossPercentage: Int,
    framesInPayloadBuffer: Int,
    packetSizeMilliseconds: Int,
) {
    val frameLength = context.frameLength
    val lpcInput = workspace.lpcInput
    if (features.signalType == SIG_TYPE_VOICED) {
        val subframeWeightsQ15 = workspace.subframeWeightsQ15
        val ltpWeights = workspace.ltpWeights
        val ltpCoefficientsQ14 = workspace.ltpCoefficientsQ14
        val correlationShifts = workspace.correlationShifts
        val ltpPredictionGainQ7 = RefInt(0)
        val residualPitchLastOffset = rightShift(frameLength, 1)
        findLtp(
            workspace.longTermPrediction,
            ltpCoefficientsQ14,
            ltpWeights,
            ltpPredictionGainQ7,
            residualPitch,
            frameLength,
            residualPitchLastOffset + frameLength,
            features.pitchL,
            subframeWeightsQ15,
            subframeLength,
            correlationShifts
        )
        features.ltPredCodGainQ7 = ltpPredictionGainQ7.value

        val periodicityIndex = RefInt(0)
        quantLtpGains(
            workspace.longTermPrediction,
            ltpCoefficientsQ14,
            features.ltpIndex,
            periodicityIndex,
            ltpWeights,
            context.ltpMuQ8,
            context.lowComplexityLtpQuantization
        )
        features.perIndex = periodicityIndex.value
        ltpCoefficientsQ14.copyInto(features.ltpCoefQ14)
        ltpScaleCtrl(
            context.longTermPredictionGains,
            packetLossPercentage,
            framesInPayloadBuffer,
            packetSizeMilliseconds,
            features,
        )
        val predictionInputOffset = frameLength - lpcOrder
        ltpAnalysisFilter(
            workspace.longTermPrediction,
            lpcInput,
            inputBuffer,
            predictionInputOffset,
            features.ltpCoefQ14,
            features.pitchL,
            inverseGainsQ16,
            subframeLength,
            lpcOrder
        )
    } else {
        features.ltpCoefQ14.fill(0)
        features.ltPredCodGainQ7 = 0

        var subframeIndex = 0
        while (subframeIndex < NB_SUBFR) {
            val inputOffset = max(
                0,
                frameLength - lpcOrder + subframeIndex * subframeLength
            )
            val outputOffset = subframeIndex * (subframeLength + lpcOrder)
            var sampleIndex = 0
            while (sampleIndex < subframeLength + lpcOrder) {
                val sample = inputBuffer[inputOffset + sampleIndex]
                lpcInput[outputOffset + sampleIndex] = smulwb(
                    inverseGainsQ16[subframeIndex],
                    sample
                )
                sampleIndex++
            }
            subframeIndex++
        }
    }
}

private fun analyzeNlsf(
    context: PredictionAnalysisContext,
    features: FrameFeatures,
    workspace: PredictionAnalysisWorkspace,
    lpcInput: IntArray,
    previousNlsfQ15: IntArray,
    nlsfQ15: IntArray,
    firstFrameAfterReset: Boolean,
    speechActivityQ8: Int,
    lpcOrder: Int,
    subframeLength: Int,
): RefInt {
    val interpolationIndex = RefInt(4)
    findLpc(
        nlsfQ15,
        interpolationIndex,
        previousNlsfQ15,
        workspace.linearPrediction,
        workspace.nlsfConversion,
        if (context.interpolatedNlsfs && !firstFrameAfterReset) 1 else 0,
        lpcOrder,
        lpcInput,
        0,
        subframeLength + lpcOrder
    )
    features.nlsfInterpCoefQ2 = interpolationIndex.value

    val nlsfWeightsQ6 = workspace.nlsfWeightsQ6
    computeLaroiaNlsfWeights(nlsfWeightsQ6, nlsfQ15, lpcOrder)
    val shouldInterpolate = context.interpolatedNlsfs && interpolationIndex.value < 4
    if (shouldInterpolate) {
        val interpolatedNlsfQ15 = workspace.interpolatedNlsfQ15
        val interpolatedWeightsQ6 = workspace.interpolatedWeightsQ6
        interpolateValues(
            interpolatedNlsfQ15,
            previousNlsfQ15,
            nlsfQ15,
            interpolationIndex.value,
            lpcOrder
        )
        computeLaroiaNlsfWeights(interpolatedWeightsQ6, interpolatedNlsfQ15, lpcOrder)

        val squaredInterpolationQ15 =
            smulbb(interpolationIndex.value, interpolationIndex.value) shl 11
        var coefficientIndex = 0
        while (coefficientIndex < lpcOrder) {
            nlsfWeightsQ6[coefficientIndex] = smlawb(
                nlsfWeightsQ6[coefficientIndex] shr 1,
                interpolatedWeightsQ6[coefficientIndex],
                squaredInterpolationQ15
            )
            coefficientIndex++
        }
    }
    val nlsfIndices = workspace.nlsfIndices
    val nlsfWeightQ15 = if (features.signalType == SIG_TYPE_VOICED) {
        smlawb(66, -8388, speechActivityQ8)
    } else {
        smlawb(164, -33554, speechActivityQ8)
    }
    val nlsfFluctuationWeightQ16 = if (features.signalType == SIG_TYPE_VOICED) {
        smlawb(6554, -838848, speechActivityQ8)
    } else {
        smlawb(13107, -1677696, speechActivityQ8 + features.sparsenessQ8)
    }

    encodeNlsfMsvq(
        workspace.nlsfQuantization,
        nlsfIndices,
        nlsfQ15,
        context.nlsfCodebooks[features.signalType],
        previousNlsfQ15,
        nlsfWeightsQ6,
        max(1, nlsfWeightQ15),
        nlsfFluctuationWeightQ16,
        context.nlsfSurvivorCount,
        lpcOrder,
        if (firstFrameAfterReset) 1 else 0
    )
    var stageIndex = 0
    while (stageIndex < nlsfIndices.size) {
        features.nlsfIndices[stageIndex] = nlsfIndices[stageIndex]
        stageIndex++
    }
    return interpolationIndex
}

private fun buildPredictionCoefficients(
    features: FrameFeatures,
    workspace: PredictionAnalysisWorkspace,
    previousNlsfQ15: IntArray,
    nlsfQ15: IntArray,
    interpolationIndex: RefInt,
    lpcOrder: Int
) {
    val firstPredictionCoefficientsQ12 = workspace.firstPredictionCoefficientsQ12
    val secondPredictionCoefficientsQ12 = workspace.secondPredictionCoefficientsQ12
    convertNlsfToStablePrediction(
        secondPredictionCoefficientsQ12,
        nlsfQ15,
        lpcOrder,
        workspace.nlsfConversion
    )
    if (interpolationIndex.value < 4) {
        val interpolatedNlsfQ15 = workspace.interpolatedNlsfQ15
        interpolateValues(
            interpolatedNlsfQ15,
            previousNlsfQ15,
            nlsfQ15,
            interpolationIndex.value,
            lpcOrder
        )
        convertNlsfToStablePrediction(
            firstPredictionCoefficientsQ12,
            interpolatedNlsfQ15,
            lpcOrder,
            workspace.nlsfConversion
        )
    } else {
        secondPredictionCoefficientsQ12.copyInto(firstPredictionCoefficientsQ12)
    }
    firstPredictionCoefficientsQ12.copyInto(features.predCoefQ12, 0, 0, lpcOrder)
    secondPredictionCoefficientsQ12.copyInto(features.predCoefQ12, MAX_LPC_ORDER, 0, lpcOrder)
}

private fun restoreZeroGains(features: FrameFeatures, subframeLength: Int) {
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        if (features.gainsQ16[subframeIndex] <= 0) {
            val energy = max(1, features.residualEnergy[subframeIndex])
            val rootMeanSquare = sqrtApprox(max(energy / max(subframeLength, 1), 1))
            features.gainsQ16[subframeIndex] = max(rootMeanSquare shl 8, 1 shl 16)
        }
        subframeIndex++
    }
}

// Ported from tsilk.
// Source: residual_energy_FIX.ts

@Suppress("SameParameterValue")
private fun residualEnergy(
    workspace: PredictionAnalysisWorkspace,
    energies: IntArray,
    energyShifts: IntArray,
    input: IntArray,
    inputOffset: Int,
    gains: IntArray,
    subframeLength: Int,
    lpcOrder: Int
) {
    val subframeStride = lpcOrder + subframeLength
    val residual = workspace.predictionResidual
    var analysisInputOffset = inputOffset

    var coefficientSetIndex = 0
    while (coefficientSetIndex < 2) {
        val filterState = workspace.predictionFilterState
        filterState.fill(0)
        lpcAnalysisFilter(
            input,
            analysisInputOffset,
            workspace.predictionPairsQ12[coefficientSetIndex],
            filterState,
            residual,
            (NB_SUBFR shr 1) * subframeStride,
            lpcOrder
        )
        var residualOffset = lpcOrder
        var halfSubframeIndex = 0
        while (halfSubframeIndex < NB_SUBFR shr 1) {
            val energyOutput = RefInt(0)
            val energyShiftOutput = RefInt(0)
            sumSquaresShift(
                energyOutput,
                energyShiftOutput,
                residual,
                residualOffset,
                subframeLength
            )

            val subframeIndex =
                coefficientSetIndex * (NB_SUBFR shr 1) + halfSubframeIndex
            energies[subframeIndex] = energyOutput.value
            energyShifts[subframeIndex] = -energyShiftOutput.value
            residualOffset += subframeStride
            halfSubframeIndex++
        }

        analysisInputOffset += (NB_SUBFR shr 1) * subframeStride
        coefficientSetIndex++
    }

    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val energyLeadingZeros = clz32(energies[subframeIndex]) - 1
        val gainLeadingZeros = clz32(gains[subframeIndex]) - 1
        var normalizedGain = lshift32(gains[subframeIndex], gainLeadingZeros)
        normalizedGain = smmul(normalizedGain, normalizedGain)
        energies[subframeIndex] = smmul(
            normalizedGain,
            lshift32(energies[subframeIndex], energyLeadingZeros)
        )
        energyShifts[subframeIndex] += energyLeadingZeros + 2 * gainLeadingZeros - 64
        subframeIndex++
    }
}

internal fun residualEnergy16Covar(
    coefficients: IntArray,
    coefficientOffset: Int,
    covarianceMatrix: IntArray,
    covarianceOffset: Int,
    crossCorrelation: IntArray,
    signalEnergy: Int,
    order: Int,
    coefficientQ: Int
): Int {
    var normalizationShifts = 16 - coefficientQ
    var extraShifts = normalizationShifts
    var maximumCoefficient = 0

    var index = 0
    while (index < order) {
        maximumCoefficient = max32(
            maximumCoefficient,
            abs(coefficients[coefficientOffset + index])
        )
        index++
    }

    extraShifts = minInt(extraShifts, clz32(maximumCoefficient) - 17)
    val maximumCovariance = max32(
        covarianceMatrix[covarianceOffset],
        covarianceMatrix[covarianceOffset + order * order - 1]
    )
    extraShifts = minInt(
        extraShifts,
        clz32(
            multiply(
                order,
                rightShift(smulwb(maximumCovariance, maximumCoefficient), 4)
            )
        ) - 5
    )
    extraShifts = maxInt(extraShifts, 0)

    val normalizedCoefficients = IntArray(order)
    index = 0
    while (index < order) {
        normalizedCoefficients[index] =
            coefficients[coefficientOffset + index] shl extraShifts
        index++
    }
    normalizationShifts -= extraShifts

    var crossTerm = 0
    index = 0
    while (index < order) {
        crossTerm = smlawb(
            crossTerm,
            crossCorrelation[index],
            normalizedCoefficients[index]
        )
        index++
    }

    var energy = rightShift(signalEnergy, 1 + normalizationShifts) - crossTerm
    var quadraticTerm = 0
    index = 0
    while (index < order) {
        var rowEnergy = 0
        val rowOffset = covarianceOffset + index * order
        var column = index + 1
        while (column < order) {
            rowEnergy = smlawb(
                rowEnergy,
                covarianceMatrix[rowOffset + column],
                normalizedCoefficients[column]
            )
            column++
        }
        rowEnergy = smlawb(
            rowEnergy,
            rightShift(covarianceMatrix[rowOffset + index], 1),
            normalizedCoefficients[index]
        )
        quadraticTerm = smlawb(
            quadraticTerm,
            rowEnergy,
            normalizedCoefficients[index]
        )
        index++
    }

    energy = addLshift32(energy, quadraticTerm, normalizationShifts)
    energy = if (energy < 1) {
        1
    } else if (energy > rightShift(Int.MAX_VALUE, normalizationShifts + 2)) {
        Int.MAX_VALUE shr 1
    } else {
        leftShift(energy, normalizationShifts + 1)
    }

    return energy
}
