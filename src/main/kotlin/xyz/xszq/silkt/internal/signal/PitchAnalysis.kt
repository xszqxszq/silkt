package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.tables.PitchTables.FLAT_CONTOUR_BIAS_Q20
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MAX_LAG_MS
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MIN_LAG_MS
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_CBKS_STAGE_2
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_CBKS_STAGE_2_EXT
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_CBKS_STAGE_3_MAX
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_STAGE_3_LAGS
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_SUB_FR
import xyz.xszq.silkt.internal.tables.PitchTables.PREVIOUS_LAG_BIAS_Q15
import xyz.xszq.silkt.internal.tables.PitchTables.SHORT_LAG_BIAS_Q15
import xyz.xszq.silkt.internal.tables.PitchTables.cbLagsStage2
import xyz.xszq.silkt.internal.tables.PitchTables.cbLagsStage3
import xyz.xszq.silkt.internal.tables.PitchTables.cbkOffsetsStage3
import xyz.xszq.silkt.internal.tables.PitchTables.cbkSizesStage3
import xyz.xszq.silkt.internal.tables.PitchTables.lagRangeStage3

internal const val MAX_FIND_PITCH_LPC_ORDER: Int = 16

private const val FIND_PITCH_WHITE_NOISE_FRACTION: Int = 66
private const val FIND_PITCH_BANDWIDTH_EXPANSION: Int = 64881

// Ported from tsilk SILK sources.

// Ported from tsilk.
// Source: pitch_analysis_core.ts

private fun maxAbsoluteValue(signal: IntArray, offset: Int, length: Int): Int {
    var max = 0
    var index = offset
    val end = offset + length
    while (index < end) {
        val value = abs(signal[index])
        if (value > max) {
            max = value
        }
        index++
    }
    return max
}

internal fun innerProductAligned(
    firstSignal: IntArray,
    firstOffset: Int,
    secondSignal: IntArray,
    secondOffset: Int,
    length: Int
): Int {
    var sum = 0
    var index = 0
    while (index < length) {
        sum += firstSignal[firstOffset + index] * secondSignal[secondOffset + index]
        index++
    }
    return sum
}

@Suppress("DuplicatedCode", "SameParameterValue")
private fun insertionSortDecreasingInt16(
    values: IntArray,
    valueOffset: Int,
    indices: IntArray,
    indexOffset: Int,
    valueCount: Int,
    topCount: Int
) {
    var index = 0
    while (index < topCount) {
        indices[indexOffset + index] = index
        index++
    }

    index = 1
    while (index < topCount) {
        val value = values[valueOffset + index]
        var insertionIndex = index - 1
        while (insertionIndex >= 0 && value > values[valueOffset + insertionIndex]) {
            values[valueOffset + insertionIndex + 1] = values[valueOffset + insertionIndex]
            indices[indexOffset + insertionIndex + 1] = indices[indexOffset + insertionIndex]
            insertionIndex--
        }
        values[valueOffset + insertionIndex + 1] = value
        indices[indexOffset + insertionIndex + 1] = index
        index++
    }

    index = topCount
    while (index < valueCount) {
        val value = values[valueOffset + index]
        if (value > values[valueOffset + topCount - 1]) {
            var insertionIndex = topCount - 2
            while (insertionIndex >= 0 && value > values[valueOffset + insertionIndex]) {
                values[valueOffset + insertionIndex + 1] = values[valueOffset + insertionIndex]
                indices[indexOffset + insertionIndex + 1] = indices[indexOffset + insertionIndex]
                insertionIndex--
            }
            values[valueOffset + insertionIndex + 1] = value
            indices[indexOffset + insertionIndex + 1] = index
        }
        index++
    }
}

@Suppress("SameParameterValue")
private fun findScaling(signal: IntArray, offset: Int, length: Int, sumSquareLength: Int): Int {
    val max = maxAbsoluteValue(signal, offset, length)
    val bits = if (max < 32767) {
        32 - clz32(smulbb(max, max))
    } else {
        30
    }
    val requiredBits = bits + 17 - clz16(sumSquareLength)
    return if (requiredBits < 31) 0 else requiredBits - 30
}

@Suppress("DuplicatedCode", "SameParameterValue")
private fun calculateStage3Correlations(
    stage3Correlations: IntArray,
    stage3LagValues: IntArray,
    signal: IntArray,
    signalOffset: Int,
    startLag: Int,
    sfLength: Int,
    complexity: Int
) {
    val codebookOffset = cbkOffsetsStage3[complexity]
    val codebookSize = cbkSizesStage3[complexity]
    var targetOffset = signalOffset + (sfLength shl 2)

    var subframe = 0
    while (subframe < PITCH_EST_NB_SUB_FR) {
        val lagRangeStart = lagRangeStage3[complexity][subframe][0]
        val lagRangeEnd = lagRangeStage3[complexity][subframe][1]

        var lagValueCount = 0
        var lag = lagRangeStart
        while (lag <= lagRangeEnd) {
            val basisOffset = targetOffset - (startLag + lag)
            stage3LagValues[lagValueCount] = innerProductAligned(signal, targetOffset, signal, basisOffset, sfLength)
            lagValueCount++
            lag++
        }

        var codebookIndex = codebookOffset
        while (codebookIndex < codebookOffset + codebookSize) {
            val lagOffset = cbLagsStage3[subframe][codebookIndex] - lagRangeStart
            var stage3Lag = 0
            while (stage3Lag < PITCH_EST_NB_STAGE_3_LAGS) {
                val outputIndex =
                    subframe * PITCH_EST_NB_CBKS_STAGE_3_MAX * PITCH_EST_NB_STAGE_3_LAGS +
                            codebookIndex * PITCH_EST_NB_STAGE_3_LAGS +
                            stage3Lag
                stage3Correlations[outputIndex] = stage3LagValues[lagOffset + stage3Lag]
                stage3Lag++
            }
            codebookIndex++
        }
        targetOffset += sfLength
        subframe++
    }
}

@Suppress("DuplicatedCode", "SameParameterValue")
private fun calculateStage3Energies(
    stage3Energies: IntArray,
    stage3LagValues: IntArray,
    signal: IntArray,
    signalOffset: Int,
    startLag: Int,
    sfLength: Int,
    complexity: Int
) {
    val codebookOffset = cbkOffsetsStage3[complexity]
    val codebookSize = cbkSizesStage3[complexity]
    var targetOffset = signalOffset + (sfLength shl 2)

    var subframe = 0
    while (subframe < PITCH_EST_NB_SUB_FR) {
        val lagRangeStart = lagRangeStage3[complexity][subframe][0]
        val lagRangeEnd = lagRangeStage3[complexity][subframe][1]
        val basisOffset = targetOffset - (startLag + lagRangeStart)
        var energy = innerProductAligned(signal, basisOffset, signal, basisOffset, sfLength)
        var lagEnergyCount = 0
        stage3LagValues[lagEnergyCount] = energy
        lagEnergyCount++

        var lagOffset = 1
        while (lagOffset < lagRangeEnd - lagRangeStart + 1) {
            energy -= smulbb(signal[basisOffset + sfLength - lagOffset], signal[basisOffset + sfLength - lagOffset])
            energy += smulbb(signal[basisOffset - lagOffset], signal[basisOffset - lagOffset])
            stage3LagValues[lagEnergyCount] = energy
            lagEnergyCount++
            lagOffset++
        }

        var codebookIndex = codebookOffset
        while (codebookIndex < codebookOffset + codebookSize) {
            val codebookLagOffset = cbLagsStage3[subframe][codebookIndex] - lagRangeStart
            var stage3Lag = 0
            while (stage3Lag < PITCH_EST_NB_STAGE_3_LAGS) {
                val outputIndex =
                    subframe * PITCH_EST_NB_CBKS_STAGE_3_MAX * PITCH_EST_NB_STAGE_3_LAGS +
                            codebookIndex * PITCH_EST_NB_STAGE_3_LAGS +
                            stage3Lag
                stage3Energies[outputIndex] = stage3LagValues[codebookLagOffset + stage3Lag]
                stage3Lag++
            }
            codebookIndex++
        }
        targetOffset += sfLength
        subframe++
    }
}

@Suppress("SameParameterValue")
private fun preparePitchSignals(
    signal: IntArray,
    workspace: PitchAnalysisWorkspace,
    sampleRateKHz: Int,
    frameLength: Int,
    frameLength8KHz: Int
) {
    val signal8KHz = workspace.signal8KHz
    val signal4KHz = workspace.signal4KHz
    val decimationState = workspace.decimationState

    when (sampleRateKHz) {
        16 -> {
            decimationState.fill(0, 0, 2)
            resamplerDown2(decimationState, 0, signal8KHz, 0, signal, 0, frameLength)
        }

        12 -> {
            decimationState.fill(0, 0, 6)
            resamplerDown23(decimationState, 0, signal8KHz, 0, signal, 0, frameLength)
        }

        24 -> {
            decimationState.fill(0, 0, 8)
            resamplerDown3(decimationState, 0, signal8KHz, 0, signal, 0, frameLength)
        }

        else -> signal.copyInto(signal8KHz, 0, 0, frameLength8KHz)
    }

    decimationState.fill(0)
    resamplerDown2(decimationState, 0, signal4KHz, 0, signal8KHz, 0, frameLength8KHz)

    var index = frameLength / 4 - 1
    while (index > 0) {
        signal4KHz[index] = addSat16(signal4KHz[index], signal4KHz[index - 1])
        index--
    }

    val sumSquareLength = max(frameLength8KHz / 8, frameLength / 4 / 2)
    val scalingShift = findScaling(signal4KHz, 0, frameLength / 4, sumSquareLength)
    if (scalingShift > 0) {
        index = 0
        while (index < frameLength / 4) {
            signal4KHz[index] = rightShift(signal4KHz[index], scalingShift)
            index++
        }
    }
}

@Suppress("SameParameterValue")
private fun calculateCorrelations4KHz(
    workspace: PitchAnalysisWorkspace,
    frameLength4KHz: Int,
    minLag4KHz: Int,
    maxLag4KHz: Int,
    subframeLength8KHz: Int
) {
    val signal4KHz = workspace.signal4KHz
    val correlations = workspace.correlations
    var targetOffset = frameLength4KHz / 2

    var subframe = 0
    while (subframe < 2) {
        var basisOffset = targetOffset - minLag4KHz
        var crossCorrelation =
            innerProductAligned(signal4KHz, targetOffset, signal4KHz, basisOffset, subframeLength8KHz)
        var normalizer = innerProductAligned(signal4KHz, basisOffset, signal4KHz, basisOffset, subframeLength8KHz)
        normalizer += subframeLength8KHz * 4000
        correlations[subframe][minLag4KHz] = sat16(div32(crossCorrelation, sqrtApprox(normalizer) + 1))

        var lag = minLag4KHz + 1
        while (lag <= maxLag4KHz) {
            basisOffset--
            crossCorrelation = innerProductAligned(
                signal4KHz,
                targetOffset,
                signal4KHz,
                basisOffset,
                subframeLength8KHz
            )
            normalizer += signal4KHz[basisOffset] * signal4KHz[basisOffset] -
                    signal4KHz[basisOffset + subframeLength8KHz] * signal4KHz[basisOffset + subframeLength8KHz]
            correlations[subframe][lag] = sat16(div32(crossCorrelation, sqrtApprox(normalizer) + 1))
            lag++
        }
        targetOffset += subframeLength8KHz
        subframe++
    }

    var lag = maxLag4KHz
    while (lag >= minLag4KHz) {
        var combinedCorrelation = correlations[0][lag] + correlations[1][lag]
        combinedCorrelation = rightShift(combinedCorrelation, 1)
        combinedCorrelation = smlawb(combinedCorrelation, combinedCorrelation, leftShift(-lag, 4))
        correlations[0][lag] = sat16(combinedCorrelation)
        lag--
    }
}

@Suppress("SameParameterValue")
private fun expandLagCandidates(
    searchLags: IntArray,
    lagCandidates: IntArray,
    seedLagCount: Int,
    minLag8KHz: Int,
    maxLag8KHz: Int
): Int {
    var index = minLag8KHz - 5
    while (index < maxLag8KHz + 5) {
        lagCandidates[index] = 0
        index++
    }
    index = 0
    while (index < seedLagCount) {
        lagCandidates[searchLags[index]] = 1
        index++
    }

    index = maxLag8KHz + 3
    while (index >= minLag8KHz) {
        lagCandidates[index] += lagCandidates[index - 1] + lagCandidates[index - 2]
        index--
    }

    var searchLagCount = 0
    index = minLag8KHz
    while (index < maxLag8KHz + 1) {
        if (lagCandidates[index + 1] > 0) {
            searchLags[searchLagCount] = index
            searchLagCount++
        }
        index++
    }

    index = maxLag8KHz + 3
    while (index >= minLag8KHz) {
        lagCandidates[index] +=
            lagCandidates[index - 1] + lagCandidates[index - 2] + lagCandidates[index - 3]
        index--
    }

    var lagCandidateCount = 0
    index = minLag8KHz
    while (index < maxLag8KHz + 4) {
        if (lagCandidates[index] > 0) {
            lagCandidates[lagCandidateCount] = index - 2
            lagCandidateCount++
        }
        index++
    }
    return lagCandidateCount
}

@Suppress("SameParameterValue")
private fun calculateCorrelations8KHz(
    workspace: PitchAnalysisWorkspace,
    frameLength4KHz: Int,
    lagCandidateCount: Int,
    subframeLength8KHz: Int
) {
    val signal8KHz = workspace.signal8KHz
    val lagCandidates = workspace.lagCandidates
    val correlations = workspace.correlations
    correlations.forEach { it.fill(0) }

    var targetOffset = frameLength4KHz
    var subframe = 0
    while (subframe < PITCH_EST_NB_SUB_FR) {
        val targetEnergy = innerProductAligned(signal8KHz, targetOffset, signal8KHz, targetOffset, subframeLength8KHz)
        var candidateIndex = 0
        while (candidateIndex < lagCandidateCount) {
            val lag = lagCandidates[candidateIndex]
            val basisOffset = targetOffset - lag
            val crossCorrelation =
                innerProductAligned(signal8KHz, targetOffset, signal8KHz, basisOffset, subframeLength8KHz)
            val basisEnergy = innerProductAligned(signal8KHz, basisOffset, signal8KHz, basisOffset, subframeLength8KHz)
            if (crossCorrelation > 0) {
                val largerEnergy = max(targetEnergy, basisEnergy)
                var leadingZeros = clz32(crossCorrelation)
                var normalizationShift = limit32(leadingZeros - 1, 0, 15)
                var normalizedCorrelation = div32(
                    leftShift(crossCorrelation, normalizationShift),
                    rightShift(largerEnergy, 15 - normalizationShift) + 1
                )
                normalizedCorrelation = smulwb(crossCorrelation, normalizedCorrelation)
                normalizedCorrelation *= 2
                leadingZeros = clz32(normalizedCorrelation)
                normalizationShift = limit32(leadingZeros - 1, 0, 15)
                val smallerEnergy = min(targetEnergy, basisEnergy)
                correlations[subframe][lag] = div32(
                    leftShift(normalizedCorrelation, normalizationShift),
                    rightShift(smallerEnergy, 15 - normalizationShift) + 1
                )
            } else {
                correlations[subframe][lag] = 0
            }
            candidateIndex++
        }
        targetOffset += subframeLength8KHz
        subframe++
    }
}

@Suppress("SameParameterValue")
private fun selectStage2Lag(
    workspace: PitchAnalysisWorkspace,
    searchLagCount: Int,
    previousLag: Int,
    searchThresholdQ15: Int,
    sampleRateKHz: Int,
    minLag8KHz: Int,
    complexity: Int,
    ltpCorrelationQ15: Int,
    forLjc: Int
): Stage2Selection? {
    val correlations = workspace.correlations
    val searchLags = workspace.searchLags
    val stage2Scores = workspace.stage2Scores
    var bestLag = -1
    var bestCodebookIndex = 0
    var bestCorrelation = Int.MIN_VALUE
    var bestBiasedScore = Int.MIN_VALUE

    var previousLagLog2Q7 = 0
    var normalizedPreviousLag = previousLag
    if (normalizedPreviousLag > 0) {
        normalizedPreviousLag = when (sampleRateKHz) {
            12 -> div3216(leftShift(normalizedPreviousLag, 1), 3)
            16 -> rightShift(normalizedPreviousLag, 1)
            24 -> div3216(normalizedPreviousLag, 3)
            else -> normalizedPreviousLag
        }
        previousLagLog2Q7 = lin2Log(normalizedPreviousLag)
    }

    val stage2CodebookCount =
        if (sampleRateKHz == 8 && complexity > 0) PITCH_EST_NB_CBKS_STAGE_2_EXT else PITCH_EST_NB_CBKS_STAGE_2

    var searchIndex = 0
    while (searchIndex < searchLagCount) {
        val lag = searchLags[searchIndex]
        var codebookIndex = 0
        while (codebookIndex < stage2CodebookCount) {
            stage2Scores[codebookIndex] = 0
            var subframe = 0
            while (subframe < PITCH_EST_NB_SUB_FR) {
                stage2Scores[codebookIndex] +=
                    correlations[subframe][lag + cbLagsStage2[subframe][codebookIndex]]
                subframe++
            }
            codebookIndex++
        }

        var candidateCorrelation = Int.MIN_VALUE
        var candidateCodebookIndex = 0
        codebookIndex = 0
        while (codebookIndex < stage2CodebookCount) {
            if (stage2Scores[codebookIndex] > candidateCorrelation) {
                candidateCorrelation = stage2Scores[codebookIndex]
                candidateCodebookIndex = codebookIndex
            }
            codebookIndex++
        }

        val lagLog2Q7 = lin2Log(lag)
        val candidateScore = if (forLjc != 0) {
            candidateCorrelation
        } else {
            candidateCorrelation -
                    rightShift(smulbb(PITCH_EST_NB_SUB_FR * SHORT_LAG_BIAS_Q15, lagLog2Q7), 7)
        }
        val biasedScore = if (previousLag > 0) {
            var deltaLagLog2SquaredQ7 = lagLog2Q7 - previousLagLog2Q7
            deltaLagLog2SquaredQ7 = rightShift(
                smulbb(deltaLagLog2SquaredQ7, deltaLagLog2SquaredQ7),
                7
            )
            var previousLagBiasQ15 = rightShift(
                smulbb(PITCH_EST_NB_SUB_FR * PREVIOUS_LAG_BIAS_Q15, ltpCorrelationQ15),
                15
            )
            previousLagBiasQ15 = div32(
                multiply(previousLagBiasQ15, deltaLagLog2SquaredQ7),
                deltaLagLog2SquaredQ7 + (1 shl 6)
            )
            candidateScore - previousLagBiasQ15
        } else {
            candidateScore
        }

        if (
            biasedScore > bestBiasedScore &&
            candidateCorrelation > searchThresholdQ15 &&
            cbLagsStage2[0][candidateCodebookIndex] <= minLag8KHz
        ) {
            bestBiasedScore = biasedScore
            bestCorrelation = candidateCorrelation
            bestLag = lag
            bestCodebookIndex = candidateCodebookIndex
        }
        searchIndex++
    }

    return if (bestLag == -1) {
        null
    } else {
        Stage2Selection(bestLag, bestCodebookIndex, bestCorrelation)
    }
}

private fun refineHighRatePitch(
    signal: IntArray,
    workspace: PitchAnalysisWorkspace,
    pitchOut: IntArray,
    pitchOutOffset: Int,
    stage2: Stage2Selection,
    sampleRateKHz: Int,
    complexity: Int,
    frameLength: Int,
    subframeLength: Int,
    minLag: Int,
    maxLag: Int
): PitchSelection {
    val scalingShift = findScaling(signal, 0, frameLength, subframeLength)
    val analysisSignal = if (scalingShift > 0) {
        val scaledSignal = workspace.analysisSignal
        var index = 0
        while (index < frameLength) {
            scaledSignal[index] = rightShift(signal[index], scalingShift)
            index++
        }
        scaledSignal
    } else {
        signal
    }

    val stage2CodebookIndex = stage2.codebookIndex
    var lag = when (sampleRateKHz) {
        12 -> rightShift(smulbb(stage2.lag, 3), 1)
        16 -> leftShift(stage2.lag, 1)
        else -> smulbb(stage2.lag, 3)
    }
    lag = limitInt(lag, minLag, maxLag)
    val startLag = max(lag - 2, minLag)
    val endLag = min(lag + 2, maxLag)
    val stage2LtpCorrelationQ15 =
        if (stage2.correlationQ15 > 0) sqrtApprox(leftShift(stage2.correlationQ15, 13)) else 0

    var subframe = 0
    while (subframe < PITCH_EST_NB_SUB_FR) {
        pitchOut[pitchOutOffset + subframe] =
            lag + 2 * cbLagsStage2[subframe][stage2CodebookIndex]
        subframe++
    }

    val stage3Correlations = workspace.stage3Correlations
    val stage3Energies = workspace.stage3Energies
    calculateStage3Correlations(
        stage3Correlations,
        workspace.stage3LagValues,
        analysisSignal,
        0,
        startLag,
        subframeLength,
        complexity
    )
    calculateStage3Energies(
        stage3Energies,
        workspace.stage3LagValues,
        analysisSignal,
        0,
        startLag,
        subframeLength,
        complexity
    )

    val contourBias = div3216(FLAT_CONTOUR_BIAS_Q20, lag)
    val codebookCount = cbkSizesStage3[complexity]
    val codebookOffset = cbkOffsetsStage3[complexity]
    var bestLag = lag
    var bestCodebookIndex = 0
    var bestScore = Int.MIN_VALUE
    var lagOffset = 0
    var candidateLag = startLag
    while (candidateLag <= endLag) {
        var codebookIndex = codebookOffset
        while (codebookIndex < codebookOffset + codebookCount) {
            var crossCorrelation = 0
            var energy = 0
            subframe = 0
            while (subframe < PITCH_EST_NB_SUB_FR) {
                val index =
                    subframe * PITCH_EST_NB_CBKS_STAGE_3_MAX * PITCH_EST_NB_STAGE_3_LAGS +
                            codebookIndex * PITCH_EST_NB_STAGE_3_LAGS +
                            lagOffset
                energy += rightShift(stage3Energies[index], 2)
                crossCorrelation += rightShift(stage3Correlations[index], 2)
                subframe++
            }

            var score = 0
            if (crossCorrelation > 0) {
                val leadingZeros = clz32(crossCorrelation)
                val normalizationShift = limit32(leadingZeros - 1, 0, 13)
                score = div32(
                    leftShift(crossCorrelation, normalizationShift),
                    rightShift(energy, 13 - normalizationShift) + 1
                )
                score = sat16(score)
                score = smulwb(crossCorrelation, score)
                score = if (score > rightShift(Int.MAX_VALUE, 3)) {
                    Int.MAX_VALUE
                } else {
                    leftShift(score, 3)
                }
                var contourDistance = codebookIndex - rightShift(PITCH_EST_NB_CBKS_STAGE_3_MAX, 1)
                contourDistance = multiply(contourDistance, contourDistance)
                val contourWeight = 32767 - rightShift(multiply(contourBias, contourDistance), 5)
                score = leftShift(smulwb(score, contourWeight), 1)
            }

            if (score > bestScore && candidateLag + cbLagsStage3[0][codebookIndex] <= maxLag) {
                bestScore = score
                bestLag = candidateLag
                bestCodebookIndex = codebookIndex
            }
            codebookIndex++
        }
        lagOffset++
        candidateLag++
    }
    return PitchSelection(bestLag, bestCodebookIndex, stage2LtpCorrelationQ15)
}

@Suppress("SameReturnValue")
private fun clearPitchAnalysisResult(
    pitchOut: IntArray,
    pitchOutOffset: Int,
    lagAndContourIndices: IntArray,
    ltpCorrelationQ15: IntArray
): Int {
    pitchOut.fill(0, pitchOutOffset, pitchOutOffset + PITCH_EST_NB_SUB_FR)
    ltpCorrelationQ15[0] = 0
    lagAndContourIndices[0] = 0
    lagAndContourIndices[1] = 0
    return 1
}

// Ported from tsilk.
// Source: pitch_analysis_core_main.ts

@Suppress("SameParameterValue")
private fun pitchAnalysisCore(
    signal: IntArray,
    workspace: PitchAnalysisWorkspace,
    pitchOut: IntArray,
    pitchOutOffset: Int,
    lagAndContourIndices: IntArray,
    ltpCorrelationQ15: IntArray,
    prevLag: Int,
    primarySearchThresholdQ16: Int,
    secondarySearchThresholdQ15: Int,
    sampleRateKHz: Int,
    complexity: Int,
    forLjc: Int
): Int {
    val signal4KHz = workspace.signal4KHz
    val correlations = workspace.correlations
    val searchLags = workspace.searchLags
    val frameLength = PITCH_EST_FRAME_LENGTH_MS * sampleRateKHz
    val frameLength4KHz = PITCH_EST_FRAME_LENGTH_MS * 4
    val frameLength8KHz = PITCH_EST_FRAME_LENGTH_MS * 8
    val subframeLength = frameLength shr 3
    val subframeLength8KHz = frameLength8KHz shr 3
    val minLag = PITCH_EST_MIN_LAG_MS * sampleRateKHz
    val minLag4KHz = PITCH_EST_MIN_LAG_MS * 4
    val minLag8KHz = PITCH_EST_MIN_LAG_MS * 8
    val maxLag = PITCH_EST_MAX_LAG_MS * sampleRateKHz
    val maxLag4KHz = PITCH_EST_MAX_LAG_MS * 4
    val maxLag8KHz = PITCH_EST_MAX_LAG_MS * 8

    preparePitchSignals(signal, workspace, sampleRateKHz, frameLength, frameLength8KHz)
    calculateCorrelations4KHz(
        workspace,
        frameLength4KHz,
        minLag4KHz,
        maxLag4KHz,
        subframeLength8KHz
    )
    var seedLagCount = 4 + 2 * complexity
    insertionSortDecreasingInt16(
        correlations[0],
        minLag4KHz,
        searchLags,
        0,
        maxLag4KHz - minLag4KHz + 1,
        seedLagCount
    )

    val targetOffset = frameLength4KHz shr 1
    var signalEnergy = innerProductAligned(signal4KHz, targetOffset, signal4KHz, targetOffset, frameLength4KHz shr 1)
    signalEnergy = addPosSat32(signalEnergy, 1000)
    val maximumCorrelation = correlations[0][minLag4KHz]
    if (rightShift(signalEnergy, 4 + 2) > maximumCorrelation * maximumCorrelation) {
        return clearPitchAnalysisResult(
            pitchOut,
            pitchOutOffset,
            lagAndContourIndices,
            ltpCorrelationQ15
        )
    }

    val correlationThreshold = smulwb(primarySearchThresholdQ16, maximumCorrelation)
    var index = 0
    while (index < seedLagCount) {
        if (correlations[0][minLag4KHz + index] > correlationThreshold) {
            searchLags[index] = (searchLags[index] + minLag4KHz) shl 1
        } else {
            seedLagCount = index
            break
        }
        index++
    }
    val lagCandidateCount = expandLagCandidates(
        searchLags,
        workspace.lagCandidates,
        seedLagCount,
        minLag8KHz,
        maxLag8KHz
    )
    val scalingShift = findScaling(workspace.signal8KHz, 0, frameLength8KHz, subframeLength8KHz)
    if (scalingShift > 0) {
        var index = 0
        while (index < frameLength8KHz) {
            workspace.signal8KHz[index] = rightShift(workspace.signal8KHz[index], scalingShift)
            index++
        }
    }
    calculateCorrelations8KHz(
        workspace,
        frameLength4KHz,
        lagCandidateCount,
        subframeLength8KHz
    )
    val stage2ThresholdQ15 = rightShift(
        smulbb(secondarySearchThresholdQ15, secondarySearchThresholdQ15),
        13
    )
    val stage2 = selectStage2Lag(
        workspace,
        seedLagCount,
        prevLag,
        stage2ThresholdQ15,
        sampleRateKHz,
        minLag8KHz,
        complexity,
        ltpCorrelationQ15[0],
        forLjc
    )
    if (stage2 == null) {
        return clearPitchAnalysisResult(
            pitchOut,
            pitchOutOffset,
            lagAndContourIndices,
            ltpCorrelationQ15
        )
    }

    val selection = if (sampleRateKHz > 8) {
        refineHighRatePitch(
            signal,
            workspace,
            pitchOut,
            pitchOutOffset,
            stage2,
            sampleRateKHz,
            complexity,
            frameLength,
            subframeLength,
            minLag,
            maxLag
        )
    } else {
        val ltpCorrelation =
            if (stage2.correlationQ15 > 0) sqrtApprox(leftShift(stage2.correlationQ15, 13)) else 0
        PitchSelection(stage2.lag, stage2.codebookIndex, ltpCorrelation)
    }

    ltpCorrelationQ15[0] = selection.ltpCorrelationQ15
    var subframe = 0
    while (subframe < PITCH_EST_NB_SUB_FR) {
        val codebookLags = if (sampleRateKHz > 8) cbLagsStage3 else cbLagsStage2
        pitchOut[pitchOutOffset + subframe] =
            selection.lag + codebookLags[subframe][selection.codebookIndex]
        subframe++
    }
    lagAndContourIndices[0] =
        selection.lag - if (sampleRateKHz > 8) minLag else minLag8KHz
    lagAndContourIndices[1] = selection.codebookIndex
    return 0
}

private fun buildWeightedPitchSignal(
    input: IntArray,
    workspace: PitchAnalysisWorkspace,
    context: PitchAnalysisContext,
    inputOffset: Int
) {
    val bufferLength = context.lookahead + (context.frameLength shl 1)
    val inputBufferOffset = inputOffset - context.frameLength
    val weightedSignal = workspace.weightedSignal
    var inputReadOffset = inputBufferOffset + bufferLength - context.windowLength
    var weightedSignalOffset = 0

    applySineWindow(weightedSignal, weightedSignalOffset, input, inputReadOffset, 1, context.lookahead)
    weightedSignalOffset += context.lookahead
    inputReadOffset += context.lookahead

    val middleLength = context.windowLength - (context.lookahead shl 1)
    var index = 0
    while (index < middleLength) {
        weightedSignal[weightedSignalOffset + index] = input[inputReadOffset + index]
        index++
    }
    weightedSignalOffset += middleLength
    inputReadOffset += middleLength
    applySineWindow(weightedSignal, weightedSignalOffset, input, inputReadOffset, 2, context.lookahead)
}

private fun analyzePitchResidual(
    input: IntArray,
    residualPitch: IntArray,
    workspace: PitchAnalysisWorkspace,
    context: PitchAnalysisContext,
    features: FrameFeatures,
    inputOffset: Int
) {
    val bufferLength = context.lookahead + (context.frameLength shl 1)
    val inputBufferOffset = inputOffset - context.frameLength
    val weightedSignal = workspace.weightedSignal
    val autoCorrelation = workspace.autoCorrelation
    val reflectionCoefficientsQ15 = workspace.reflectionCoefficientsQ15
    val coefficientsQ24 = workspace.coefficientsQ24
    val predictionState = workspace.predictionState
    val coefficientsQ12 = workspace.coefficientsQ12

    autocorr(
        autoCorrelation,
        workspace.autocorrelationScale,
        weightedSignal,
        0,
        context.windowLength,
        context.lpcOrder + 1
    )
    autoCorrelation[0] =
        smlawb(autoCorrelation[0], autoCorrelation[0], FIND_PITCH_WHITE_NOISE_FRACTION)
    val residualEnergy = schur(reflectionCoefficientsQ15, autoCorrelation, context.lpcOrder)
    features.predGainQ16 = div32VarQ(autoCorrelation[0], max(residualEnergy, 1), 16)
    k2A(coefficientsQ24, reflectionCoefficientsQ15, context.lpcOrder)

    var index = 0
    while (index < context.lpcOrder) {
        coefficientsQ12[index] = sat16(rightShift(coefficientsQ24[index], 12))
        index++
    }
    expandBandwidth(coefficientsQ12, 0, context.lpcOrder, FIND_PITCH_BANDWIDTH_EXPANSION)
    predictionState.fill(0)
    movingAveragePrediction(
        input,
        inputBufferOffset,
        coefficientsQ12,
        0,
        predictionState,
        0,
        residualPitch,
        0,
        bufferLength,
        context.lpcOrder
    )
    residualPitch.fill(0, 0, context.lpcOrder)
}

private fun calculatePitchSearchThresholdQ15(
    context: PitchAnalysisContext,
    features: FrameFeatures,
    speechActivityQ8: Int
): Int {
    var searchThresholdQ15 = 14745
    searchThresholdQ15 = smlabb(searchThresholdQ15, -131, context.lpcOrder)
    searchThresholdQ15 = smlabb(searchThresholdQ15, -12, speechActivityQ8)
    searchThresholdQ15 = smlabb(searchThresholdQ15, 4915, context.previousSignalType)
    searchThresholdQ15 = smlawb(searchThresholdQ15, -6554, features.inputTiltQ15)
    return sat16(searchThresholdQ15)
}

// Ported from tsilk.
// Source: find_pitch_lags_FIX.ts

internal fun findPitchLags(
    context: PitchAnalysisContext,
    features: FrameFeatures,
    speechActivityQ8: Int,
    residualPitch: IntArray,
    input: IntArray,
    inputOffset: Int
) {
    val workspace = context.workspace
    buildWeightedPitchSignal(input, workspace, context, inputOffset)
    analyzePitchResidual(input, residualPitch, workspace, context, features, inputOffset)
    val searchThresholdQ15 =
        calculatePitchSearchThresholdQ15(context, features, speechActivityQ8)
    val lagAndContourIndices = workspace.lagAndContourIndices
    val ltpCorrelation = workspace.ltpCorrelationOutput
    lagAndContourIndices[0] = features.lagIndex
    lagAndContourIndices[1] = features.contourIndex
    ltpCorrelation[0] = context.ltpCorrelationQ15
    features.signalType = pitchAnalysisCore(
        residualPitch,
        workspace,
        features.pitchL,
        0,
        lagAndContourIndices,
        ltpCorrelation,
        context.previousLag,
        context.correlationThresholdQ16,
        searchThresholdQ15,
        context.sampleRateKHz,
        context.complexity,
        0
    )
    features.lagIndex = lagAndContourIndices[0]
    features.contourIndex = lagAndContourIndices[1]
    context.ltpCorrelationQ15 = ltpCorrelation[0]
}
