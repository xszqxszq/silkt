package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER

// Ported from tsilk.
// Source: nlsf_vq.ts

private const val MIN_NLSF_SPACING = 3
private const val WEIGHT_FRACTION_BITS = 6
private const val WEIGHT_SCALE = 1 shl (15 + WEIGHT_FRACTION_BITS)
private const val MAX_WEIGHT_Q6 = 32767

private const val MAX_MSVQ_SURVIVORS = 16
private const val MAX_MSVQ_STAGES = 10
private const val MAX_RATE_DISTORTIONS = 256
private const val PRUNING_FACTOR_Q16 = 6554

internal class NlsfQuantizationWorkspace(
    val rateDistortionQ18: IntArray,
    val accumulatedRatesQ5: IntArray,
    val nextAccumulatedRatesQ5: IntArray,
    val sortedCandidateIndexes: IntArray,
    val survivorPaths: IntArray,
    val nextSurvivorPaths: IntArray,
    val residualQ15: IntArray,
    val nextResidualQ15: IntArray,
    val candidatePath: IntArray,
    val packedWeightsQ6: IntArray,
) {
    constructor() : this(
        rateDistortionQ18 = IntArray(MAX_RATE_DISTORTIONS),
        accumulatedRatesQ5 = IntArray(MAX_MSVQ_SURVIVORS),
        nextAccumulatedRatesQ5 = IntArray(MAX_MSVQ_SURVIVORS),
        sortedCandidateIndexes = IntArray(MAX_MSVQ_SURVIVORS),
        survivorPaths = IntArray(MAX_MSVQ_SURVIVORS * MAX_MSVQ_STAGES),
        nextSurvivorPaths = IntArray(MAX_MSVQ_SURVIVORS * MAX_MSVQ_STAGES),
        residualQ15 = IntArray(MAX_MSVQ_SURVIVORS * MAX_MSVQ_SURVIVORS),
        nextResidualQ15 = IntArray(MAX_MSVQ_SURVIVORS * MAX_MSVQ_SURVIVORS),
        candidatePath = IntArray(MAX_MSVQ_STAGES),
        packedWeightsQ6 = IntArray(MAX_LPC_ORDER shr 1),
    )
}

internal fun computeLaroiaNlsfWeights(
    weightsQ6: IntArray,
    nlsfQ15: IntArray,
    order: Int
) {
    var leftSpacingWeight = WEIGHT_SCALE / maxInt(nlsfQ15[0], MIN_NLSF_SPACING)
    var rightSpacingWeight =
        WEIGHT_SCALE / maxInt(nlsfQ15[1] - nlsfQ15[0], MIN_NLSF_SPACING)
    weightsQ6[0] = minInt(leftSpacingWeight + rightSpacingWeight, MAX_WEIGHT_Q6)

    var coefficientIndex = 1
    while (coefficientIndex < order - 1) {
        leftSpacingWeight = rightSpacingWeight
        rightSpacingWeight = WEIGHT_SCALE / maxInt(
            nlsfQ15[coefficientIndex + 1] - nlsfQ15[coefficientIndex],
            MIN_NLSF_SPACING
        )
        weightsQ6[coefficientIndex] =
            minInt(leftSpacingWeight + rightSpacingWeight, MAX_WEIGHT_Q6)

        leftSpacingWeight = rightSpacingWeight
        rightSpacingWeight = WEIGHT_SCALE / maxInt(
            nlsfQ15[coefficientIndex + 2] - nlsfQ15[coefficientIndex + 1],
            MIN_NLSF_SPACING
        )
        weightsQ6[coefficientIndex + 1] =
            minInt(leftSpacingWeight + rightSpacingWeight, MAX_WEIGHT_Q6)
        coefficientIndex += 2
    }

    leftSpacingWeight = rightSpacingWeight
    rightSpacingWeight = WEIGHT_SCALE / maxInt(
        (1 shl 15) - nlsfQ15[order - 1],
        MIN_NLSF_SPACING
    )
    weightsQ6[order - 1] = minInt(rightSpacingWeight + leftSpacingWeight, MAX_WEIGHT_Q6)
}

internal fun interpolateValues(
    interpolated: IntArray,
    start: IntArray,
    end: IntArray,
    interpolationFactorQ2: Int,
    length: Int
) {
    var valueIndex = 0
    while (valueIndex < length) {
        interpolated[valueIndex] = start[valueIndex] + rightShift(
            smulbb(
                end[valueIndex] - start[valueIndex],
                interpolationFactorQ2
            ),
            2
        )
        valueIndex++
    }
}

private fun calculateNlsfVectorErrors(
    errorsQ20: IntArray,
    packedWeightsQ6: IntArray,
    inputVectorsQ15: IntArray,
    inputOffset: Int,
    weightsQ6: IntArray,
    codebookVectorsQ15: IntArray,
    inputVectorCount: Int,
    codebookVectorCount: Int,
    order: Int
) {
    if (order == 16) {
        calculateNlsfVectorErrorsOrder16(
            errorsQ20,
            packedWeightsQ6,
            inputVectorsQ15,
            inputOffset,
            codebookVectorsQ15,
            inputVectorCount,
            codebookVectorCount
        )
        return
    }
    if (order == 10) {
        calculateNlsfVectorErrorsOrder10(
            errorsQ20,
            packedWeightsQ6,
            inputVectorsQ15,
            inputOffset,
            codebookVectorsQ15,
            inputVectorCount,
            codebookVectorCount
        )
        return
    }

    var errorOffset = 0
    var vectorInputOffset = inputOffset
    var inputVectorIndex = 0
    while (inputVectorIndex < inputVectorCount) {
        var codebookOffset = 0
        var codebookVectorIndex = 0
        while (codebookVectorIndex < codebookVectorCount) {
            var weightedErrorQ20 = 0
            var coefficientIndex = 0
            while (coefficientIndex < order) {
                val packedWeightQ6 = packedWeightsQ6[coefficientIndex shr 1]
                var differenceQ15 =
                    inputVectorsQ15[vectorInputOffset + coefficientIndex] -
                            codebookVectorsQ15[codebookOffset++]
                weightedErrorQ20 = smlawb(
                    weightedErrorQ20,
                    smulbb(differenceQ15, differenceQ15),
                    packedWeightQ6
                )
                differenceQ15 =
                    inputVectorsQ15[vectorInputOffset + coefficientIndex + 1] -
                            codebookVectorsQ15[codebookOffset++]
                weightedErrorQ20 = smlawt(
                    weightedErrorQ20,
                    smulbb(differenceQ15, differenceQ15),
                    packedWeightQ6
                )
                coefficientIndex += 2
            }
            errorsQ20[errorOffset + codebookVectorIndex] = weightedErrorQ20
            codebookVectorIndex++
        }
        errorOffset += codebookVectorCount
        vectorInputOffset += order
        inputVectorIndex++
    }
}

@Suppress("DuplicatedCode")
private fun calculateNlsfVectorErrorsOrder16(
    errorsQ20: IntArray,
    packedWeightsQ6: IntArray,
    inputVectorsQ15: IntArray,
    inputOffset: Int,
    codebookVectorsQ15: IntArray,
    inputVectorCount: Int,
    codebookVectorCount: Int
) {
    var errorOffset = 0
    var vectorInputOffset = inputOffset
    var inputVectorIndex = 0
    while (inputVectorIndex < inputVectorCount) {
        var codebookOffset = 0
        var codebookVectorIndex = 0
        while (codebookVectorIndex < codebookVectorCount) {
            var weightedErrorQ20 = 0
            var pairIndex = 0
            while (pairIndex < 8) {
                val packedWeightQ6 = packedWeightsQ6[pairIndex]
                var differenceQ15 =
                    inputVectorsQ15[vectorInputOffset] - codebookVectorsQ15[codebookOffset]
                weightedErrorQ20 = smlawb(
                    weightedErrorQ20,
                    smulbb(differenceQ15, differenceQ15),
                    packedWeightQ6
                )
                differenceQ15 =
                    inputVectorsQ15[vectorInputOffset + 1] - codebookVectorsQ15[codebookOffset + 1]
                weightedErrorQ20 = smlawt(
                    weightedErrorQ20,
                    smulbb(differenceQ15, differenceQ15),
                    packedWeightQ6
                )
                vectorInputOffset += 2
                codebookOffset += 2
                pairIndex++
            }
            errorsQ20[errorOffset + codebookVectorIndex] = weightedErrorQ20
            codebookVectorIndex++
            vectorInputOffset -= 16
        }
        errorOffset += codebookVectorCount
        vectorInputOffset += 16
        inputVectorIndex++
    }
}

@Suppress("DuplicatedCode")
private fun calculateNlsfVectorErrorsOrder10(
    errorsQ20: IntArray,
    packedWeightsQ6: IntArray,
    inputVectorsQ15: IntArray,
    inputOffset: Int,
    codebookVectorsQ15: IntArray,
    inputVectorCount: Int,
    codebookVectorCount: Int
) {
    var errorOffset = 0
    var vectorInputOffset = inputOffset
    var inputVectorIndex = 0
    while (inputVectorIndex < inputVectorCount) {
        var codebookOffset = 0
        var codebookVectorIndex = 0
        while (codebookVectorIndex < codebookVectorCount) {
            var weightedErrorQ20 = 0
            var pairIndex = 0
            while (pairIndex < 5) {
                val packedWeightQ6 = packedWeightsQ6[pairIndex]
                var differenceQ15 =
                    inputVectorsQ15[vectorInputOffset] - codebookVectorsQ15[codebookOffset]
                weightedErrorQ20 = smlawb(
                    weightedErrorQ20,
                    smulbb(differenceQ15, differenceQ15),
                    packedWeightQ6
                )
                differenceQ15 =
                    inputVectorsQ15[vectorInputOffset + 1] - codebookVectorsQ15[codebookOffset + 1]
                weightedErrorQ20 = smlawt(
                    weightedErrorQ20,
                    smulbb(differenceQ15, differenceQ15),
                    packedWeightQ6
                )
                vectorInputOffset += 2
                codebookOffset += 2
                pairIndex++
            }
            errorsQ20[errorOffset + codebookVectorIndex] = weightedErrorQ20
            codebookVectorIndex++
            vectorInputOffset -= 10
        }
        errorOffset += codebookVectorCount
        vectorInputOffset += 10
        inputVectorIndex++
    }
}

private fun packNlsfWeights(
    packedWeightsQ6: IntArray,
    weightsQ6: IntArray,
    order: Int
) {
    var packedWeightIndex = 0
    while (packedWeightIndex < order shr 1) {
        packedWeightsQ6[packedWeightIndex] = weightsQ6[2 * packedWeightIndex] or
                (weightsQ6[2 * packedWeightIndex + 1] shl 16)
        packedWeightIndex++
    }
}

@Suppress("SameParameterValue")
private fun calculateNlsfRateDistortion(
    rateDistortionQ20: IntArray,
    packedWeightsQ6: IntArray,
    codebookVectorsQ15: IntArray,
    codebookRatesQ5: IntArray,
    codebookVectorCount: Int,
    inputVectorsQ15: IntArray,
    inputOffset: Int,
    weightsQ6: IntArray,
    accumulatedRatesQ5: IntArray,
    rateWeightQ15: Int,
    inputVectorCount: Int,
    order: Int
) {
    calculateNlsfVectorErrors(
        rateDistortionQ20,
        packedWeightsQ6,
        inputVectorsQ15,
        inputOffset,
        weightsQ6,
        codebookVectorsQ15,
        inputVectorCount,
        codebookVectorCount,
        order
    )

    var rateDistortionOffset = 0
    var inputVectorIndex = 0
    while (inputVectorIndex < inputVectorCount) {
        var codebookVectorIndex = 0
        while (codebookVectorIndex < codebookVectorCount) {
            rateDistortionQ20[rateDistortionOffset + codebookVectorIndex] = smlabb(
                rateDistortionQ20[rateDistortionOffset + codebookVectorIndex],
                accumulatedRatesQ5[inputVectorIndex] +
                        codebookRatesQ5[codebookVectorIndex],
                rateWeightQ15
            )
            codebookVectorIndex++
        }
        rateDistortionOffset += codebookVectorCount
        inputVectorIndex++
    }
}

@Suppress("DuplicatedCode")
private fun sortSmallestWithIndexes(
    values: IntArray,
    indexes: IntArray,
    valueCount: Int,
    survivorCount: Int
) {
    var index = 0
    while (index < survivorCount) {
        indexes[index] = index
        index++
    }

    index = 1
    while (index < survivorCount) {
        val value = values[index]
        val sourceIndex = indexes[index]
        var insertionIndex = index - 1
        while (insertionIndex >= 0 && value < values[insertionIndex]) {
            values[insertionIndex + 1] = values[insertionIndex]
            indexes[insertionIndex + 1] = indexes[insertionIndex]
            insertionIndex--
        }
        values[insertionIndex + 1] = value
        indexes[insertionIndex + 1] = sourceIndex
        index++
    }

    index = survivorCount
    while (index < valueCount) {
        val value = values[index]
        if (value < values[survivorCount - 1]) {
            var insertionIndex = survivorCount - 2
            while (insertionIndex >= 0 && value < values[insertionIndex]) {
                values[insertionIndex + 1] = values[insertionIndex]
                indexes[insertionIndex + 1] = indexes[insertionIndex]
                insertionIndex--
            }
            values[insertionIndex + 1] = value
            indexes[insertionIndex + 1] = index
        }
        index++
    }
}

@Suppress("DuplicatedCode")
internal fun encodeNlsfMsvq(
    workspace: NlsfQuantizationWorkspace,
    nlsfIndices: IntArray,
    nlsfQ15: IntArray,
    codebook: NlsfCodebook,
    previousQuantizedNlsfQ15: IntArray,
    weightsQ6: IntArray,
    nlsfRateWeightQ15: Int,
    fluctuationReductionWeightQ16: Int,
    maximumSurvivorCount: Int,
    order: Int,
    fluctuationReductionDisabled: Int
) {
    val minimumSurvivorCount = maximumSurvivorCount / 2
    val rateDistortionQ18 = workspace.rateDistortionQ18
    val accumulatedRatesQ5 = workspace.accumulatedRatesQ5
    val nextAccumulatedRatesQ5 = workspace.nextAccumulatedRatesQ5
    val sortedCandidateIndexes = workspace.sortedCandidateIndexes
    val survivorPaths = workspace.survivorPaths
    val nextSurvivorPaths = workspace.nextSurvivorPaths
    val residualQ15 = workspace.residualQ15
    val nextResidualQ15 = workspace.nextResidualQ15
    packNlsfWeights(workspace.packedWeightsQ6, weightsQ6, order)
    accumulatedRatesQ5[0] = 0
    nlsfQ15.copyInto(residualQ15, 0, 0, order)

    var previousSurvivorCount = 1
    var stageIndex = 0
    while (stageIndex < codebook.stageCount) {
        val stage = codebook.stages[stageIndex]
        var currentSurvivorCount = min(
            maximumSurvivorCount,
            previousSurvivorCount * stage.vectorCount
        )
        calculateNlsfRateDistortion(
            rateDistortionQ18,
            workspace.packedWeightsQ6,
            stage.vectorsQ15,
            stage.ratesQ5,
            stage.vectorCount,
            residualQ15,
            0,
            weightsQ6,
            accumulatedRatesQ5,
            nlsfRateWeightQ15,
            previousSurvivorCount,
            order
        )
        sortSmallestWithIndexes(
            rateDistortionQ18,
            sortedCandidateIndexes,
            previousSurvivorCount * stage.vectorCount,
            currentSurvivorCount
        )

        currentSurvivorCount = pruneMsvqSurvivors(
            rateDistortionQ18,
            currentSurvivorCount,
            minimumSurvivorCount,
            maximumSurvivorCount
        )

        var survivorIndex = 0
        while (survivorIndex < currentSurvivorCount) {
            var inputSurvivorIndex = 0
            var codebookVectorIndex: Int
            if (stageIndex > 0) {
                val candidateIndex = sortedCandidateIndexes[survivorIndex]
                if (stage.vectorCount == 8) {
                    inputSurvivorIndex = candidateIndex shr 3
                    codebookVectorIndex = candidateIndex and 7
                } else {
                    inputSurvivorIndex = candidateIndex / stage.vectorCount
                    codebookVectorIndex = candidateIndex % stage.vectorCount
                }
            } else {
                codebookVectorIndex = sortedCandidateIndexes[survivorIndex]
            }

            val inputVectorOffset = inputSurvivorIndex * order
            val codebookVectorOffset = codebookVectorIndex * order
            val survivorResidualOffset = survivorIndex * order
            var coefficientIndex = 0
            while (coefficientIndex < order) {
                nextResidualQ15[survivorResidualOffset + coefficientIndex] =
                    residualQ15[inputVectorOffset + coefficientIndex] -
                            stage.vectorsQ15[codebookVectorOffset + coefficientIndex]
                coefficientIndex++
            }
            nextAccumulatedRatesQ5[survivorIndex] =
                accumulatedRatesQ5[inputSurvivorIndex] +
                        stage.ratesQ5[codebookVectorIndex]

            val inputPathOffset = inputSurvivorIndex * codebook.stageCount
            val survivorPathOffset = survivorIndex * codebook.stageCount
            var previousStageIndex = 0
            while (previousStageIndex < stageIndex) {
                nextSurvivorPaths[survivorPathOffset + previousStageIndex] =
                    survivorPaths[inputPathOffset + previousStageIndex]
                previousStageIndex++
            }
            nextSurvivorPaths[survivorPathOffset + stageIndex] = codebookVectorIndex
            survivorIndex++
        }

        if (stageIndex < codebook.stageCount - 1) {
            nextResidualQ15.copyInto(residualQ15)
            nextAccumulatedRatesQ5.copyInto(accumulatedRatesQ5)
            nextSurvivorPaths.copyInto(survivorPaths)
        }
        previousSurvivorCount = currentSurvivorCount
        stageIndex++
    }

    val bestSurvivorIndex =
        if (fluctuationReductionDisabled == 1) {
            0
        } else {
            selectFluctuationReducedSurvivor(
                nlsfQ15,
                previousQuantizedNlsfQ15,
                weightsQ6,
                rateDistortionQ18,
                nextSurvivorPaths,
                workspace.candidatePath,
                codebook,
                fluctuationReductionWeightQ16,
                previousSurvivorCount,
                order
            )
        }

    val bestPathOffset = bestSurvivorIndex * codebook.stageCount
    nextSurvivorPaths.copyInto(
        nlsfIndices,
        0,
        bestPathOffset,
        bestPathOffset + codebook.stageCount
    )
    decodeNlsfMsvq(nlsfQ15, codebook, nlsfIndices, order)
}

private fun pruneMsvqSurvivors(
    rateDistortionQ18: IntArray,
    currentSurvivorCount: Int,
    minimumSurvivorCount: Int,
    maximumSurvivorCount: Int
): Int {
    var survivorCount = currentSurvivorCount
    if (rateDistortionQ18[0] < Int.MAX_VALUE / 16) {
        val distortionLimit = smlawb(
            rateDistortionQ18[0],
            maximumSurvivorCount * rateDistortionQ18[0],
            PRUNING_FACTOR_Q16
        )
        while (
            rateDistortionQ18[survivorCount - 1] > distortionLimit &&
            survivorCount > minimumSurvivorCount
        ) {
            survivorCount--
        }
    }
    return survivorCount
}

private fun selectFluctuationReducedSurvivor(
    nlsfQ15: IntArray,
    previousQuantizedNlsfQ15: IntArray,
    weightsQ6: IntArray,
    rateDistortionQ18: IntArray,
    survivorPaths: IntArray,
    candidatePath: IntArray,
    codebook: NlsfCodebook,
    fluctuationReductionWeightQ16: Int,
    survivorCount: Int,
    order: Int
): Int {
    var bestSurvivorIndex = 0
    var bestFluctuationCostQ20 = Int.MAX_VALUE
    var survivorIndex = 0
    while (survivorIndex < survivorCount) {
        val survivorPathOffset = survivorIndex * codebook.stageCount
        survivorPaths.copyInto(
            candidatePath,
            0,
            survivorPathOffset,
            survivorPathOffset + codebook.stageCount
        )
        decodeNlsfMsvq(nlsfQ15, codebook, candidatePath, order)

        var weightedSquaredErrorQ20 = 0
        var coefficientIndex = 0
        while (coefficientIndex < order) {
            var differenceQ15 =
                nlsfQ15[coefficientIndex] - previousQuantizedNlsfQ15[coefficientIndex]
            weightedSquaredErrorQ20 = smlawb(
                weightedSquaredErrorQ20,
                smulbb(differenceQ15, differenceQ15),
                weightsQ6[coefficientIndex]
            )
            differenceQ15 =
                nlsfQ15[coefficientIndex + 1] - previousQuantizedNlsfQ15[coefficientIndex + 1]
            weightedSquaredErrorQ20 = smlawb(
                weightedSquaredErrorQ20,
                smulbb(differenceQ15, differenceQ15),
                weightsQ6[coefficientIndex + 1]
            )
            coefficientIndex += 2
        }

        val fluctuationCostQ20 = addPosSat32(
            rateDistortionQ18[survivorIndex],
            smulwb(weightedSquaredErrorQ20, fluctuationReductionWeightQ16)
        )
        if (fluctuationCostQ20 < bestFluctuationCostQ20) {
            bestFluctuationCostQ20 = fluctuationCostQ20
            bestSurvivorIndex = survivorIndex
        }
        survivorIndex++
    }
    return bestSurvivorIndex
}
