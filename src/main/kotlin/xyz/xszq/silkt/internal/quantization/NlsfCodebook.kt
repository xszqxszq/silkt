package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.fixedpoint.max
import xyz.xszq.silkt.internal.fixedpoint.min
import xyz.xszq.silkt.internal.fixedpoint.rshiftRound
import xyz.xszq.silkt.internal.fixedpoint.toInt32

// Ported from tsilk.
// Source: nlsf.ts

internal class NlsfCodebook(
    val stageCount: Int,
    val stages: Array<NlsfCodebookStage>,
    val minimumDeltasQ15: IntArray,
    val cdf: IntArray,
    val cdfStartOffsets: IntArray
)

internal fun decodeNlsfMsvq(
    nlsfQ15: IntArray,
    codebook: NlsfCodebook,
    nlsfIndices: IntArray,
    lpcOrder: Int
) {
    val firstStageVectors = codebook.stages[0].vectorsQ15
    val firstVectorOffset = nlsfIndices[0] * lpcOrder
    var coefficientIndex = 0
    while (coefficientIndex < lpcOrder) {
        nlsfQ15[coefficientIndex] = firstStageVectors[firstVectorOffset + coefficientIndex]
        coefficientIndex++
    }

    var stageIndex = 1
    while (stageIndex < codebook.stageCount) {
        val stageVectors = codebook.stages[stageIndex].vectorsQ15
        val vectorOffset = if (lpcOrder == 16) {
            nlsfIndices[stageIndex] shl 4
        } else {
            nlsfIndices[stageIndex] * lpcOrder
        }
        coefficientIndex = 0
        while (coefficientIndex < lpcOrder) {
            nlsfQ15[coefficientIndex] =
                toInt32(nlsfQ15[coefficientIndex] + stageVectors[vectorOffset + coefficientIndex])
            coefficientIndex++
        }
        stageIndex++
    }
    nlsfStabilize(nlsfQ15, codebook.minimumDeltasQ15, lpcOrder)
}

private const val MAX_STABILIZATION_ATTEMPTS = 20

private fun nlsfStabilize(nlsfQ15: IntArray, minimumDeltasQ15: IntArray, lpcOrder: Int) {
    var attempt = 0
    var coefficientIndex: Int
    while (attempt < MAX_STABILIZATION_ATTEMPTS) {
        var minimumDifferenceQ15 = nlsfQ15[0] - minimumDeltasQ15[0]
        var violationIndex = 0
        coefficientIndex = 1
        while (coefficientIndex <= lpcOrder - 1) {
            val difference = nlsfQ15[coefficientIndex] -
                    (nlsfQ15[coefficientIndex - 1] + minimumDeltasQ15[coefficientIndex])
            if (difference < minimumDifferenceQ15) {
                minimumDifferenceQ15 = difference
            }
            coefficientIndex++
        }

        val finalDifference = (1 shl 15) - (nlsfQ15[lpcOrder - 1] + minimumDeltasQ15[lpcOrder])
        if (finalDifference < minimumDifferenceQ15) {
            minimumDifferenceQ15 = finalDifference
            violationIndex = lpcOrder
        }
        if (minimumDifferenceQ15 >= 0) {
            return
        }

        when (violationIndex) {
            0 -> {
                nlsfQ15[0] = minimumDeltasQ15[0]
            }

            lpcOrder -> {
                nlsfQ15[lpcOrder - 1] = (1 shl 15) - minimumDeltasQ15[lpcOrder]
            }

            else -> {
                var minimumCenter = 0
                var lowerIndex = 0
                while (lowerIndex < violationIndex) {
                    minimumCenter += minimumDeltasQ15[lowerIndex]
                    lowerIndex++
                }
                minimumCenter += minimumDeltasQ15[violationIndex] shr 1

                var maximumCenter = 1 shl 15
                var upperIndex = lpcOrder
                while (upperIndex > violationIndex) {
                    maximumCenter -= minimumDeltasQ15[upperIndex]
                    upperIndex--
                }
                maximumCenter -= minimumDeltasQ15[violationIndex] - (minimumDeltasQ15[violationIndex] shr 1)

                var center = rshiftRound(nlsfQ15[violationIndex - 1] + nlsfQ15[violationIndex], 1)
                center = max(minimumCenter, min(maximumCenter, center))
                nlsfQ15[violationIndex - 1] = center - (minimumDeltasQ15[violationIndex] shr 1)
                nlsfQ15[violationIndex] = nlsfQ15[violationIndex - 1] + minimumDeltasQ15[violationIndex]
            }
        }
        attempt++
    }

    insertionSort(nlsfQ15, lpcOrder)
    nlsfQ15[0] = max(nlsfQ15[0], minimumDeltasQ15[0])
    coefficientIndex = 1
    while (coefficientIndex < lpcOrder) {
        nlsfQ15[coefficientIndex] = max(
            nlsfQ15[coefficientIndex],
            nlsfQ15[coefficientIndex - 1] + minimumDeltasQ15[coefficientIndex]
        )
        coefficientIndex++
    }
    nlsfQ15[lpcOrder - 1] = min(nlsfQ15[lpcOrder - 1], (1 shl 15) - minimumDeltasQ15[lpcOrder])
    coefficientIndex = lpcOrder - 2
    while (coefficientIndex >= 0) {
        nlsfQ15[coefficientIndex] = min(
            nlsfQ15[coefficientIndex],
            nlsfQ15[coefficientIndex + 1] - minimumDeltasQ15[coefficientIndex + 1]
        )
        coefficientIndex--
    }
}

private fun insertionSort(values: IntArray, length: Int) {
    var insertionIndex = 1
    while (insertionIndex < length) {
        val value = values[insertionIndex]
        var comparisonIndex = insertionIndex - 1
        while (comparisonIndex >= 0 && values[comparisonIndex] > value) {
            values[comparisonIndex + 1] = values[comparisonIndex]
            comparisonIndex--
        }
        values[comparisonIndex + 1] = value
        insertionIndex++
    }
}

// Ported from tsilk.
// Source: nlsf_codebooks.ts

internal fun buildNlsfCodebook(
    stageCount: Int,
    stageVectorCounts: IntArray,
    cdf: IntArray,
    cdfStartOffsets: IntArray,
    minimumDeltasQ15: IntArray,
    vectorsQ15: IntArray,
    ratesQ5: IntArray,
    lpcOrder: Int
): NlsfCodebook {
    val stages = mutableListOf<NlsfCodebookStage>()
    var vectorOffset = 0
    var rateOffset = 0
    var stageIndex = 0
    while (stageIndex < stageCount) {
        val vectorCount = stageVectorCounts[stageIndex]
        val vectorDataLength = vectorCount * lpcOrder
        stages += NlsfCodebookStage(
            vectorCount = vectorCount,
            vectorsQ15 = vectorsQ15.copyOfRange(vectorOffset, vectorOffset + vectorDataLength),
            ratesQ5 = ratesQ5.copyOfRange(rateOffset, rateOffset + vectorCount)
        )
        vectorOffset += vectorDataLength
        rateOffset += vectorCount
        stageIndex++
    }
    return NlsfCodebook(
        stageCount = stageCount,
        stages = stages.toTypedArray(),
        minimumDeltasQ15 = minimumDeltasQ15,
        cdf = cdf,
        cdfStartOffsets = cdfStartOffsets
    )
}
