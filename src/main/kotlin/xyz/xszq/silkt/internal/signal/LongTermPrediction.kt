package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.model.LTP_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.tables.LtpTables.LTP_GAIN_MIDDLE_AVG_RD_Q14
import xyz.xszq.silkt.internal.tables.LtpTables.ltpGainBitsQ6Tables
import xyz.xszq.silkt.internal.tables.LtpTables.quantLtpVqSizes
import xyz.xszq.silkt.internal.tables.LtpTables.quantLtpVqTablesQ14
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.ltpScalesTableQ14
import xyz.xszq.silkt.internal.util.RefInt

// Ported from tsilk.
// Source: find_Ltp_FIX.ts

private const val LTP_CORRELATION_HEADROOM = 2

private const val LTP_DAMPING_DIV_3_Q16 = 218

@Suppress("SameParameterValue")
private fun fitLtp(
    sourceCoefficientsQ16: IntArray,
    destinationCoefficientsQ14: IntArray,
    sourceOffset: Int,
    destinationOffset: Int
) {
    var coefficientIndex = 0
    while (coefficientIndex < LTP_ORDER) {
        destinationCoefficientsQ14[destinationOffset + coefficientIndex] = sat16(
            rshiftRound(sourceCoefficientsQ16[sourceOffset + coefficientIndex], 2)
        )
        coefficientIndex++
    }
}

@Suppress("SameParameterValue")
private fun scaleVector32Q26Shift18(data: IntArray, dataOffset: Int, gainQ26: Int, dataLength: Int) {
    var sampleIndex = 0
    while (sampleIndex < dataLength) {
        val sampleOffset = dataOffset + sampleIndex
        data[sampleOffset] = ((data[sampleOffset].toLong() * gainQ26.toLong()) shr 8).toInt()
        sampleIndex++
    }
}

private fun maximumCorrelationShift(correlationShifts: IntArray): Int {
    var maximumShifts = 0
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        maximumShifts = maxInt(correlationShifts[subframeIndex], maximumShifts)
        subframeIndex++
    }
    return maximumShifts
}

private fun calculateLtpCodingGainQ7(
    signalEnergies: IntArray,
    filteredResidualEnergies: IntArray,
    correlationShifts: IntArray,
    subframeWeightsQ15: IntArray,
    maximumShifts: Int
): Int {
    var predictionResidualEnergy = 0
    var ltpResidualEnergy = 0
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        predictionResidualEnergy = add32(
            predictionResidualEnergy,
            rightShift(
                add32(
                    smulwb(signalEnergies[subframeIndex], subframeWeightsQ15[subframeIndex]),
                    1
                ),
                1 + maximumShifts - correlationShifts[subframeIndex]
            )
        )
        ltpResidualEnergy = add32(
            ltpResidualEnergy,
            rightShift(
                add32(
                    smulwb(
                        filteredResidualEnergies[subframeIndex],
                        subframeWeightsQ15[subframeIndex]
                    ),
                    1
                ),
                1 + maximumShifts - correlationShifts[subframeIndex]
            )
        )
        subframeIndex++
    }
    ltpResidualEnergy = max(ltpResidualEnergy, 1)
    val energyRatioQ16 = div32VarQ(predictionResidualEnergy, ltpResidualEnergy, 16)
    return smulbb(3, lin2Log(energyRatioQ16) - (16 shl 7))
}

private fun fillCoefficientSums(
    ltpCoefficientsQ14: IntArray,
    coefficientSumsQ14: IntArray
) {
    var ltpCoefficientOffset = 0
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        coefficientSumsQ14[subframeIndex] = 0
        var coefficientIndex = 0
        while (coefficientIndex < LTP_ORDER) {
            coefficientSumsQ14[subframeIndex] +=
                ltpCoefficientsQ14[ltpCoefficientOffset + coefficientIndex]
            coefficientIndex++
        }
        ltpCoefficientOffset += LTP_ORDER
        subframeIndex++
    }
}

private fun calculateMeanCoefficientDeltaQ12(
    coefficientSumsQ14: IntArray,
    diagonalWeights: IntArray,
    correlationShifts: IntArray,
    maximumShifts: Int
): Int {
    var maximumCoefficientSum = 0
    var maximumWeightBits = 0
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        maximumCoefficientSum = max32(
            maximumCoefficientSum,
            abs(coefficientSumsQ14[subframeIndex])
        )
        maximumWeightBits = max32(
            maximumWeightBits,
            32 - clz32(diagonalWeights[subframeIndex]) +
                    correlationShifts[subframeIndex] - maximumShifts
        )
        subframeIndex++
    }

    var extraShifts = maximumWeightBits + 32 - clz32(maximumCoefficientSum) - 14
    extraShifts -= 32 - 1 - 2 + maximumShifts
    extraShifts = maxInt(extraShifts, 0)
    val maximumShiftsWithExtra = maximumShifts + extraShifts
    var normalizer = rightShift(262, maximumShifts + extraShifts) + 1
    var weightedDelta = 0
    subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val relativeShift = maximumShiftsWithExtra - correlationShifts[subframeIndex]
        normalizer = add32(normalizer, rightShift(diagonalWeights[subframeIndex], relativeShift))
        weightedDelta = add32(
            weightedDelta,
            leftShift(
                smulww(
                    rightShift(diagonalWeights[subframeIndex], relativeShift),
                    coefficientSumsQ14[subframeIndex]
                ),
                2
            )
        )
        subframeIndex++
    }
    return div32VarQ(weightedDelta, normalizer, 12)
}

private fun stabilizeLtpCoefficients(
    ltpCoefficientsQ14: IntArray,
    coefficientDeltasQ14: IntArray,
    coefficientSumsQ14: IntArray,
    diagonalWeights: IntArray,
    correlationShifts: IntArray,
    maximumShifts: Int
) {
    val meanDeltaQ12 = calculateMeanCoefficientDeltaQ12(
        coefficientSumsQ14,
        diagonalWeights,
        correlationShifts,
        maximumShifts
    )

    var ltpCoefficientOffset = 0
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val shiftedWeight = if (2 - correlationShifts[subframeIndex] > 0) {
            rightShift(diagonalWeights[subframeIndex], 2 - correlationShifts[subframeIndex])
        } else {
            lshiftSat32(diagonalWeights[subframeIndex], correlationShifts[subframeIndex] - 2)
        }
        val gainQ26 = multiply(
            div32(6710886, rightShift(6710886, 10) + shiftedWeight),
            lshiftSat32(
                subSat32(
                    meanDeltaQ12,
                    rightShift(coefficientSumsQ14[subframeIndex], 2)
                ),
                4
            )
        )

        var deltaSum = 0
        var coefficientIndex = 0
        while (coefficientIndex < LTP_ORDER) {
            coefficientDeltasQ14[coefficientIndex] = max16(
                ltpCoefficientsQ14[ltpCoefficientOffset + coefficientIndex],
                1638
            )
            deltaSum += coefficientDeltasQ14[coefficientIndex]
            coefficientIndex++
        }

        val scaledDelta = div32(gainQ26, deltaSum)
        coefficientIndex = 0
        while (coefficientIndex < LTP_ORDER) {
            val coefficientOffset = ltpCoefficientOffset + coefficientIndex
            ltpCoefficientsQ14[coefficientOffset] = limit32(
                ltpCoefficientsQ14[coefficientOffset] +
                        smulwb(lshiftSat32(scaledDelta, 4), coefficientDeltasQ14[coefficientIndex]),
                -16000,
                28000
            )
            coefficientIndex++
        }
        ltpCoefficientOffset += LTP_ORDER
        subframeIndex++
    }
}

@Suppress("DuplicatedCode")
internal fun findLtp(
    workspace: LongTermPredictionWorkspace,
    ltpCoefficientsQ14: IntArray,
    ltpWeights: IntArray,
    ltPredCodGainQ7: RefInt?,
    residualPitch: IntArray,
    firstResidualPitchOffset: Int,
    lastResidualPitchOffset: Int,
    subframeLags: IntArray,
    subframeWeightsQ15: IntArray,
    subframeLength: Int,
    correlationShifts: IntArray
) {
    val fittedCoefficientsQ16 = workspace.fittedCoefficientsQ16
    val coefficientDeltasQ14 = workspace.coefficientDeltasQ14
    val coefficientSumsQ14 = workspace.coefficientSumsQ14
    val filteredResidualEnergies = workspace.filteredResidualEnergies
    val diagonalWeights = workspace.diagonalWeights
    val crossCorrelation = workspace.crossCorrelation
    val signalEnergies = workspace.signalEnergies
    val correlationShiftOutput = workspace.correlationShiftOutput

    var ltpCoefficientOffset = 0
    var correlationOffset = 0
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val residualPitchBaseOffset =
            if (subframeIndex >= NB_SUBFR shr 1) lastResidualPitchOffset else firstResidualPitchOffset
        val currentResidualOffset =
            residualPitchBaseOffset + (subframeIndex % (NB_SUBFR shr 1)) * subframeLength
        val lagOffset = currentResidualOffset - subframeLags[subframeIndex] - LTP_ORDER / 2
        val energyOutput = workspace.energyOutput
        val shiftOutput = workspace.shiftOutput

        sumSquaresShift(
            energyOutput,
            shiftOutput,
            residualPitch,
            currentResidualOffset,
            subframeLength
        )
        signalEnergies[subframeIndex] = energyOutput.value
        var residualShifts = shiftOutput.value
        val leadingZeros = clz32(signalEnergies[subframeIndex])
        if (leadingZeros < LTP_CORRELATION_HEADROOM) {
            signalEnergies[subframeIndex] = rshiftRound(
                signalEnergies[subframeIndex],
                LTP_CORRELATION_HEADROOM - leadingZeros
            )
            residualShifts += LTP_CORRELATION_HEADROOM - leadingZeros
        }

        correlationShifts[subframeIndex] = residualShifts
        correlationShiftOutput.value = correlationShifts[subframeIndex]
        corrMatrix(
            energyOutput,
            shiftOutput,
            residualPitch,
            lagOffset,
            subframeLength,
            LTP_ORDER,
            LTP_CORRELATION_HEADROOM,
            ltpWeights,
            correlationOffset,
            correlationShiftOutput
        )
        correlationShifts[subframeIndex] = correlationShiftOutput.value
        corrVector(
            residualPitch,
            lagOffset,
            residualPitch,
            currentResidualOffset,
            subframeLength,
            LTP_ORDER,
            crossCorrelation,
            correlationShifts[subframeIndex]
        )
        if (correlationShifts[subframeIndex] > residualShifts) {
            signalEnergies[subframeIndex] = rightShift(
                signalEnergies[subframeIndex],
                correlationShifts[subframeIndex] - residualShifts
            )
        }

        var regularization = 1
        regularization = smlawb(
            regularization,
            signalEnergies[subframeIndex],
            LTP_DAMPING_DIV_3_Q16
        )
        regularization = smlawb(
            regularization,
            ltpWeights[correlationOffset],
            LTP_DAMPING_DIV_3_Q16
        )
        regularization = smlawb(
            regularization,
            ltpWeights[correlationOffset + (LTP_ORDER - 1) * LTP_ORDER + LTP_ORDER - 1],
            LTP_DAMPING_DIV_3_Q16
        )

        val regularizedEnergy = workspace.regularizedEnergy
        regularizedEnergy.value = signalEnergies[subframeIndex]
        regularizeCorrelations(
            ltpWeights,
            correlationOffset,
            regularizedEnergy,
            regularization,
            LTP_ORDER
        )
        signalEnergies[subframeIndex] = regularizedEnergy.value

        solveLDL(
            workspace.cholesky,
            ltpWeights,
            correlationOffset,
            LTP_ORDER,
            crossCorrelation,
            fittedCoefficientsQ16
        )
        fitLtp(
            fittedCoefficientsQ16,
            ltpCoefficientsQ14,
            0,
            ltpCoefficientOffset
        )
        filteredResidualEnergies[subframeIndex] = residualEnergy16Covar(
            ltpCoefficientsQ14,
            ltpCoefficientOffset,
            ltpWeights,
            correlationOffset,
            crossCorrelation,
            signalEnergies[subframeIndex],
            LTP_ORDER,
            14
        )

        val normalizationShifts =
            minInt(correlationShifts[subframeIndex], LTP_CORRELATION_HEADROOM)
        var denominator = lshiftSat32(
            smulwb(filteredResidualEnergies[subframeIndex], subframeWeightsQ15[subframeIndex]),
            1 + normalizationShifts
        ) + rightShift(
            smulwb(subframeLength, 655),
            correlationShifts[subframeIndex] - normalizationShifts
        )
        denominator = max(denominator, 1)
        var scalingGain = div32(leftShift(subframeWeightsQ15[subframeIndex], 16), denominator)
        scalingGain = rightShift(
            scalingGain,
            31 + correlationShifts[subframeIndex] - normalizationShifts - 26
        )

        var maximumWeight = 0
        var weightIndex = 0
        while (weightIndex < LTP_ORDER * LTP_ORDER) {
            maximumWeight = max(ltpWeights[correlationOffset + weightIndex], maximumWeight)
            weightIndex++
        }
        val availableBits = clz32(maximumWeight) - 1 - 3
        if (26 - 18 + availableBits < 31) {
            scalingGain = min(scalingGain, leftShift(1, 26 - 18 + availableBits))
        }
        scaleVector32Q26Shift18(
            ltpWeights,
            correlationOffset,
            scalingGain,
            LTP_ORDER * LTP_ORDER
        )

        val centerWeightIndex = (LTP_ORDER shr 1) * LTP_ORDER + (LTP_ORDER shr 1)
        diagonalWeights[subframeIndex] = ltpWeights[correlationOffset + centerWeightIndex]
        ltpCoefficientOffset += LTP_ORDER
        correlationOffset += LTP_ORDER * LTP_ORDER
        subframeIndex++
    }

    val maximumShifts = maximumCorrelationShift(correlationShifts)
    if (ltPredCodGainQ7 != null) {
        ltPredCodGainQ7.value = calculateLtpCodingGainQ7(
            signalEnergies,
            filteredResidualEnergies,
            correlationShifts,
            subframeWeightsQ15,
            maximumShifts
        )
    }

    fillCoefficientSums(ltpCoefficientsQ14, coefficientSumsQ14)
    stabilizeLtpCoefficients(
        ltpCoefficientsQ14,
        coefficientDeltasQ14,
        coefficientSumsQ14,
        diagonalWeights,
        correlationShifts,
        maximumShifts
    )
}

// Ported from tsilk.
// Source: Ltp_analysis_filter_FIX.ts

internal fun ltpAnalysisFilter(
    workspace: LongTermPredictionWorkspace,
    residual: IntArray,
    input: IntArray,
    inputOffset: Int,
    ltpCoefficientsQ14: IntArray,
    pitchLags: IntArray,
    inverseGainsQ16: IntArray,
    subframeLength: Int,
    preLength: Int
) {
    var frameInputOffset = inputOffset
    var residualWriteOffset = 0
    val coefficientsQ14 = workspace.analysisCoefficientsQ14
    val centerIndex = LTP_ORDER shr 1

    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        var lagInputOffset = frameInputOffset - pitchLags[subframeIndex]

        var coefficientIndex = 0
        while (coefficientIndex < LTP_ORDER) {
            coefficientsQ14[coefficientIndex] =
                ltpCoefficientsQ14[subframeIndex * LTP_ORDER + coefficientIndex]
            coefficientIndex++
        }

        var sampleIndex = 0
        while (sampleIndex < subframeLength + preLength) {
            val inputSample = input[frameInputOffset + sampleIndex]
            residual[residualWriteOffset + sampleIndex] = inputSample

            var prediction = smulbb(input[lagInputOffset + centerIndex], coefficientsQ14[0])
            var coefficientOrder = 1
            while (coefficientOrder < LTP_ORDER) {
                prediction = smlabbOvflw(
                    prediction,
                    input[lagInputOffset + centerIndex - coefficientOrder],
                    coefficientsQ14[coefficientOrder]
                )
                coefficientOrder++
            }
            prediction = rshiftRound(prediction, 14)
            var sampleResidual = inputSample - prediction
            sampleResidual = sampleResidual.coerceIn(-32768, 32767)
            residual[residualWriteOffset + sampleIndex] = sampleResidual
            residual[residualWriteOffset + sampleIndex] = smulwb(
                inverseGainsQ16[subframeIndex],
                residual[residualWriteOffset + sampleIndex]
            )
            lagInputOffset++
            sampleIndex++
        }

        residualWriteOffset += subframeLength + preLength
        frameInputOffset += subframeLength
        subframeIndex++
    }
}

// Ported from tsilk.
// Source: Ltp_scale_ctrl_FIX.ts

private const val THRESHOLD_COUNT = 11

private val ltpScaleThresholdsQ15: IntArray =
    intArrayOf(31129, 26214, 16384, 13107, 9830, 6554, 4915, 3276, 2621, 2458, 0)

internal class LongTermPredictionGainState(
    var smoothedQ7: Int = 0,
    var previousQ7: Int = 0,
)

internal fun ltpScaleCtrl(
    gainState: LongTermPredictionGainState,
    packetLossPercentage: Int,
    payloadFrameCount: Int,
    packetSizeMs: Int,
    features: FrameFeatures,
) {
    gainState.smoothedQ7 =
        maxInt(features.ltPredCodGainQ7 - gainState.previousQ7, 0) +
                rshiftRound(gainState.smoothedQ7, 1)
    gainState.previousQ7 = features.ltPredCodGainQ7

    val gainOutputQ5 = rshiftRound(
        rightShift(features.ltPredCodGainQ7, 1) + rightShift(gainState.smoothedQ7, 1),
        3
    )
    val gainLimitQ15 = sigmoidQ15(gainOutputQ5 - (3 shl 5))
    var ltpScaleIndex = 0

    if (payloadFrameCount == 0) {
        val framesPerPacket = div3216(packetSizeMs, FRAME_LENGTH_MS)
        val roundLoss = packetLossPercentage + framesPerPacket - 1
        val primaryThresholdQ15 = ltpScaleThresholdsQ15[minInt(roundLoss, THRESHOLD_COUNT - 1)]
        val secondaryThresholdQ15 =
            ltpScaleThresholdsQ15[minInt(roundLoss + 1, THRESHOLD_COUNT - 1)]

        if (gainLimitQ15 > primaryThresholdQ15) {
            ltpScaleIndex = 2
        } else if (gainLimitQ15 > secondaryThresholdQ15) {
            ltpScaleIndex = 1
        }
    }

    features.ltpScaleIndex = ltpScaleIndex
    features.ltpScaleQ14 = ltpScalesTableQ14[ltpScaleIndex]
}

// Ported from tsilk.
// Source: quant_Ltp_gains_FIX.ts

internal fun quantLtpGains(
    workspace: LongTermPredictionWorkspace,
    coefficientsQ14: IntArray,
    codebookIndices: IntArray,
    periodicityIndex: RefInt,
    weightsQ18: IntArray,
    complexityWeightQ8: Int,
    lowComplexity: Boolean
) {
    val stageIndices = workspace.quantizationStageIndices
    val subframeRateDistortion = workspace.quantizationRateDistortionQ14
    val stageIndex = workspace.quantizationBestIndex
    var minimumRateDistortion = Int.MAX_VALUE

    var periodicityTableIndex = 0
    while (periodicityTableIndex < 3) {
        val stageBitsQ6 = ltpGainBitsQ6Tables[periodicityTableIndex]
        val codebookQ14 = quantLtpVqTablesQ14[periodicityTableIndex]
        val codebookSize = quantLtpVqSizes[periodicityTableIndex]
        var weightOffset = 0
        var coefficientOffset = 0
        var rateDistortion = 0

        var subframeIndex = 0
        while (subframeIndex < NB_SUBFR) {
            weightedMatrixVectorQuantize(
                stageIndex,
                subframeRateDistortion,
                coefficientsQ14,
                coefficientOffset,
                weightsQ18,
                weightOffset,
                codebookQ14,
                0,
                stageBitsQ6,
                complexityWeightQ8,
                codebookSize
            )
            stageIndices[subframeIndex] = stageIndex.value
            rateDistortion = addPosSat32(rateDistortion, subframeRateDistortion.value)
            coefficientOffset += LTP_ORDER
            weightOffset += LTP_ORDER * LTP_ORDER
            subframeIndex++
        }

        rateDistortion = min(Int.MAX_VALUE - 1, rateDistortion)
        if (rateDistortion < minimumRateDistortion) {
            minimumRateDistortion = rateDistortion

            var resultIndex = 0
            while (resultIndex < NB_SUBFR) {
                codebookIndices[resultIndex] = stageIndices[resultIndex]
                resultIndex++
            }
            periodicityIndex.value = periodicityTableIndex
        }

        if (lowComplexity && rateDistortion < LTP_GAIN_MIDDLE_AVG_RD_Q14) {
            break
        }
        periodicityTableIndex++
    }

    val selectedCodebookQ14 = quantLtpVqTablesQ14[periodicityIndex.value]
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        var coefficientIndex = 0
        while (coefficientIndex < LTP_ORDER) {
            coefficientsQ14[subframeIndex * LTP_ORDER + coefficientIndex] =
                selectedCodebookQ14[
                    mla(coefficientIndex, codebookIndices[subframeIndex], LTP_ORDER)
                ]
            coefficientIndex++
        }
        subframeIndex++
    }
}

// Ported from tsilk.
// Source: VQ_WMat_EC_FIX.ts

@Suppress("SameParameterValue")
private fun weightedMatrixVectorQuantize(
    bestIndex: RefInt,
    rateDistortionQ14: RefInt,
    inputQ14: IntArray,
    inputOffset: Int,
    weightsQ18: IntArray,
    weightOffset: Int,
    codebookQ14: IntArray,
    codebookOffset: Int,
    codeLengthsQ6: IntArray,
    complexityWeightQ8: Int,
    codebookSize: Int
) {
    rateDistortionQ14.value = Int.MAX_VALUE
    var codebookRowOffset = codebookOffset

    var codeIndex = 0
    while (codeIndex < codebookSize) {
        val packedDifference01 =
            ((inputQ14[inputOffset] - codebookQ14[codebookRowOffset]) and 65535) or
                    (((inputQ14[inputOffset + 1] - codebookQ14[codebookRowOffset + 1]) and 65535) shl 16)
        val packedDifference23 =
            ((inputQ14[inputOffset + 2] - codebookQ14[codebookRowOffset + 2]) and 65535) or
                    (((inputQ14[inputOffset + 3] - codebookQ14[codebookRowOffset + 3]) and 65535) shl 16)
        val difference4 = inputQ14[inputOffset + 4] - codebookQ14[codebookRowOffset + 4]

        var candidateRateDistortionQ14 = smulbb(complexityWeightQ8, codeLengthsQ6[codeIndex])
        var weightedSumQ16 = smulwt(weightsQ18[weightOffset + 1], packedDifference01)
        weightedSumQ16 = smlawb(weightedSumQ16, weightsQ18[weightOffset + 2], packedDifference23)
        weightedSumQ16 = smlawt(weightedSumQ16, weightsQ18[weightOffset + 3], packedDifference23)
        weightedSumQ16 = smlawb(weightedSumQ16, weightsQ18[weightOffset + 4], difference4)
        weightedSumQ16 = leftShift(weightedSumQ16, 1)
        weightedSumQ16 = smlawb(weightedSumQ16, weightsQ18[weightOffset], packedDifference01)
        candidateRateDistortionQ14 = smlawb(
            candidateRateDistortionQ14,
            weightedSumQ16,
            packedDifference01
        )

        weightedSumQ16 = smulwb(weightsQ18[weightOffset + 7], packedDifference23)
        weightedSumQ16 = smlawt(weightedSumQ16, weightsQ18[weightOffset + 8], packedDifference23)
        weightedSumQ16 = smlawb(weightedSumQ16, weightsQ18[weightOffset + 9], difference4)
        weightedSumQ16 = leftShift(weightedSumQ16, 1)
        weightedSumQ16 = smlawt(weightedSumQ16, weightsQ18[weightOffset + 6], packedDifference01)
        candidateRateDistortionQ14 = smlawt(
            candidateRateDistortionQ14,
            weightedSumQ16,
            packedDifference01
        )

        weightedSumQ16 = smulwt(weightsQ18[weightOffset + 13], packedDifference23)
        weightedSumQ16 = smlawb(weightedSumQ16, weightsQ18[weightOffset + 14], difference4)
        weightedSumQ16 = leftShift(weightedSumQ16, 1)
        weightedSumQ16 = smlawb(weightedSumQ16, weightsQ18[weightOffset + 12], packedDifference23)
        candidateRateDistortionQ14 = smlawb(
            candidateRateDistortionQ14,
            weightedSumQ16,
            packedDifference23
        )

        weightedSumQ16 = smulwb(weightsQ18[weightOffset + 19], difference4)
        weightedSumQ16 = leftShift(weightedSumQ16, 1)
        weightedSumQ16 = smlawt(weightedSumQ16, weightsQ18[weightOffset + 18], packedDifference23)
        candidateRateDistortionQ14 = smlawt(
            candidateRateDistortionQ14,
            weightedSumQ16,
            packedDifference23
        )

        weightedSumQ16 = smulwb(weightsQ18[weightOffset + 24], difference4)
        candidateRateDistortionQ14 = smlawb(
            candidateRateDistortionQ14,
            weightedSumQ16,
            difference4
        )

        if (candidateRateDistortionQ14 < rateDistortionQ14.value) {
            rateDistortionQ14.value = candidateRateDistortionQ14
            bestIndex.value = codeIndex
        }

        codebookRowOffset += 5
        codeIndex++
    }
}
