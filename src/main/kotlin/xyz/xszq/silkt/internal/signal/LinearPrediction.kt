package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.util.RefInt
import kotlin.math.floor
import kotlin.math.pow

private const val MAX_SCRATCH_LPC_ORDER = 24

internal class LinearPredictionInverseResult(
    var invGainQ30: Int = 0,
    var unstable: Boolean = false,
)

// Ported from tsilk SILK sources.

// Ported from tsilk.
// Source: autocorr.ts

internal fun innerProduct16(
    firstInput: IntArray,
    firstOffset: Int,
    secondInput: IntArray,
    secondOffset: Int,
    length: Int
): Double {
    var sum = 0.0
    var sampleIndex = 0
    while (sampleIndex < length) {
        sum += firstInput[firstOffset + sampleIndex] * secondInput[secondOffset + sampleIndex]
        sampleIndex++
    }
    return sum
}

private fun countLeadingZeros64(value: Double): Int {
    if (value == 0.0) {
        return 64
    }

    var remainingValue = value
    var leadingZeros = 0
    while (remainingValue > 0) {
        remainingValue = floor(remainingValue / 2)
        leadingZeros++
    }
    return 64 - leadingZeros
}

internal fun autocorr(
    results: IntArray,
    outputScale: RefInt,
    inputData: IntArray,
    inputDataOffset: Int,
    inputDataSize: Int,
    correlationCount: Int
) {
    val actualCorrelationCount = min(inputDataSize, correlationCount)
    val energy = innerProduct16(
        inputData,
        inputDataOffset,
        inputData,
        inputDataOffset,
        inputDataSize
    ) + 1
    val leadingZeros = countLeadingZeros64(energy)
    val rightShiftCount = 35 - leadingZeros
    outputScale.value = rightShiftCount

    if (rightShiftCount <= 0) {
        results[0] = (energy * 2.0.pow(-rightShiftCount)).toInt()
        var lag = 1
        while (lag < actualCorrelationCount) {
            val correlation = innerProduct16(
                inputData,
                inputDataOffset,
                inputData,
                inputDataOffset + lag,
                inputDataSize - lag
            )
            results[lag] = (correlation * 2.0.pow(-rightShiftCount)).toInt()
            lag++
        }
    } else {
        results[0] = floor(energy / 2.0.pow(rightShiftCount)).toInt()
        var lag = 1
        while (lag < actualCorrelationCount) {
            val correlation = innerProduct16(
                inputData,
                inputDataOffset,
                inputData,
                inputDataOffset + lag,
                inputDataSize - lag
            )
            results[lag] = floor(correlation / 2.0.pow(rightShiftCount)).toInt()
            lag++
        }
    }
}

// Ported from tsilk.
// Source: bwexpander.ts

@Suppress("DuplicatedCode")
internal fun expandBandwidth32(
    coefficients: IntArray,
    order: Int,
    chirpQ16: Int
) {
    var chirp = chirpQ16

    var coefficientIndex = 0
    while (coefficientIndex < order - 1) {
        coefficients[coefficientIndex] =
            toInt32(smulww(chirp, coefficients[coefficientIndex]))
        chirp = toInt32(smulww(chirp, chirpQ16))
        coefficientIndex++
    }
    coefficients[order - 1] = toInt32(smulww(chirp, coefficients[order - 1]))
}

// Ported from tsilk.
// Source: k2a.ts

internal fun k2A(
    predictionCoefficientsQ24: IntArray,
    reflectionCoefficientsQ15: IntArray,
    order: Int
) {
    val previousCoefficientsQ24 = IntArray(MAX_SCRATCH_LPC_ORDER)

    var coefficientOrder = 0
    while (coefficientOrder < order) {
        var coefficientIndex = 0
        while (coefficientIndex < coefficientOrder) {
            previousCoefficientsQ24[coefficientIndex] =
                predictionCoefficientsQ24[coefficientIndex]
            coefficientIndex++
        }

        coefficientIndex = 0
        while (coefficientIndex < coefficientOrder) {
            predictionCoefficientsQ24[coefficientIndex] = smlawb(
                predictionCoefficientsQ24[coefficientIndex],
                leftShift(
                    previousCoefficientsQ24[coefficientOrder - coefficientIndex - 1],
                    1
                ),
                reflectionCoefficientsQ15[coefficientOrder]
            )
            coefficientIndex++
        }
        predictionCoefficientsQ24[coefficientOrder] =
            -leftShift(reflectionCoefficientsQ15[coefficientOrder], 9)
        coefficientOrder++
    }
}

// Ported from tsilk.
// Source: k2a_Q16.ts

internal fun k2AQ16(
    predictionCoefficientsQ24: IntArray,
    reflectionCoefficientsQ16: IntArray,
    order: Int
) {
    val previousCoefficientsQ24 = IntArray(MAX_SCRATCH_LPC_ORDER)

    var coefficientOrder = 0
    while (coefficientOrder < order) {
        var coefficientIndex = 0
        while (coefficientIndex < coefficientOrder) {
            previousCoefficientsQ24[coefficientIndex] =
                predictionCoefficientsQ24[coefficientIndex]
            coefficientIndex++
        }

        coefficientIndex = 0
        while (coefficientIndex < coefficientOrder) {
            predictionCoefficientsQ24[coefficientIndex] = smlaww(
                predictionCoefficientsQ24[coefficientIndex],
                previousCoefficientsQ24[coefficientOrder - coefficientIndex - 1],
                reflectionCoefficientsQ16[coefficientOrder]
            )
            coefficientIndex++
        }
        predictionCoefficientsQ24[coefficientOrder] =
            -leftShift(reflectionCoefficientsQ16[coefficientOrder], 8)
        coefficientOrder++
    }
}

// Ported from tsilk.
// Source: lpc_analysis_filter.ts

internal fun lpcAnalysisFilter(
    input: IntArray,
    inputOffset: Int,
    coefficients: IntArray,
    state: IntArray,
    output: IntArray,
    length: Int,
    order: Int
) {
    val orderHalf = order shr 1
    var firstState: Int
    var secondState: Int
    var predictionQ12: Int
    var outputSample: Int

    var sampleIndex = 0
    while (sampleIndex < length) {
        val inputValue = input[inputOffset + sampleIndex]
        firstState = state[0]
        predictionQ12 = 0

        var coefficientIndex = 0
        while (coefficientIndex < orderHalf - 1) {
            val pairedCoefficientIndex = (coefficientIndex shl 1) + 1
            secondState = state[pairedCoefficientIndex]
            state[pairedCoefficientIndex] = firstState
            predictionQ12 = smlabb(
                predictionQ12,
                firstState,
                coefficients[pairedCoefficientIndex - 1]
            )
            predictionQ12 = smlabb(
                predictionQ12,
                secondState,
                coefficients[pairedCoefficientIndex]
            )
            firstState = state[pairedCoefficientIndex + 1]
            state[pairedCoefficientIndex + 1] = secondState
            coefficientIndex++
        }

        secondState = state[order - 1]
        state[order - 1] = firstState
        predictionQ12 = smlabb(
            predictionQ12,
            firstState,
            coefficients[order - 2]
        )
        predictionQ12 = smlabb(
            predictionQ12,
            secondState,
            coefficients[order - 1]
        )
        predictionQ12 = subSat32(inputValue shl 12, predictionQ12)
        outputSample = rshiftRound(predictionQ12, 12)
        output[sampleIndex] = when {
            outputSample > 32767 -> 32767
            outputSample < -32768 -> -32768
            else -> outputSample
        }
        state[0] = inputValue
        sampleIndex++
    }
}

// Ported from tsilk.
// Source: lpc_inv_pred_gain.ts

private const val LPC_QA = 16
private const val LPC_COEFFICIENT_LIMIT = 65520
private const val INT32_MAX = 2147483647

@Suppress("DuplicatedCode")
private fun lpcInversePredGainQA(
    coefficientsQA: Array<IntArray>,
    order: Int
): LinearPredictionInverseResult {
    var currentCoefficientsQA = coefficientsQA[order and 1]
    var inverseGainQ30 = 1 shl 30

    var coefficientOrder = order - 1
    while (coefficientOrder > 0) {
        val currentCoefficient = currentCoefficientsQA[coefficientOrder]
        if (currentCoefficient > LPC_COEFFICIENT_LIMIT ||
            currentCoefficient < -LPC_COEFFICIENT_LIMIT
        ) {
            return LinearPredictionInverseResult(invGainQ30 = 0, unstable = true)
        }

        val reflectionQ31 =
            toInt32(-leftShift(currentCoefficient, 31 - LPC_QA))
        val stabilityFactorQ30 = toInt32(
            (INT32_MAX shr 1) - smmul(reflectionQ31, reflectionQ31)
        )
        if (stabilityFactorQ30 <= 1 shl 15) {
            return LinearPredictionInverseResult(invGainQ30 = 0, unstable = true)
        }

        val inverseStabilityFactorQ16 = inverse32VarQ(stabilityFactorQ30, 46)
        inverseGainQ30 = toInt32(
            leftShift(smmul(inverseGainQ30, stabilityFactorQ30), 2)
        )

        val previousCoefficientsQA = currentCoefficientsQA
        currentCoefficientsQA = coefficientsQA[coefficientOrder and 1]
        val headroom = clz32(inverseStabilityFactorQ16) - 1
        val shiftedInverseFactorQ16 =
            leftShift(inverseStabilityFactorQ16, headroom)

        var coefficientIndex = 0
        while (coefficientIndex < coefficientOrder) {
            val updatedCoefficientQA = toInt32(
                previousCoefficientsQA[coefficientIndex] -
                        leftShift(
                            smmul(
                                previousCoefficientsQA[
                                    coefficientOrder - coefficientIndex - 1
                                ],
                                reflectionQ31
                            ),
                            1
                        )
            )
            currentCoefficientsQA[coefficientIndex] = toInt32(
                leftShift(
                    smmul(updatedCoefficientQA, shiftedInverseFactorQ16),
                    16 - headroom
                )
            )
            coefficientIndex++
        }
        coefficientOrder--
    }

    val currentCoefficient = currentCoefficientsQA[0]
    if (currentCoefficient > LPC_COEFFICIENT_LIMIT ||
        currentCoefficient < -LPC_COEFFICIENT_LIMIT
    ) {
        return LinearPredictionInverseResult(invGainQ30 = 0, unstable = true)
    }

    val reflectionQ31 = toInt32(-leftShift(currentCoefficient, 31 - LPC_QA))
    val stabilityFactorQ30 = toInt32(
        (INT32_MAX shr 1) - smmul(reflectionQ31, reflectionQ31)
    )
    inverseGainQ30 = toInt32(
        leftShift(smmul(inverseGainQ30, stabilityFactorQ30), 2)
    )
    return LinearPredictionInverseResult(
        invGainQ30 = inverseGainQ30,
        unstable = false
    )
}

internal fun lpcInversePredGain(
    predictionCoefficientsQ12: IntArray,
    order: Int
): LinearPredictionInverseResult {
    val workCoefficientsQA = arrayOf(IntArray(16), IntArray(16))
    val currentCoefficientsQA = workCoefficientsQA[order and 1]

    var coefficientIndex = 0
    while (coefficientIndex < order) {
        currentCoefficientsQA[coefficientIndex] =
            leftShift(predictionCoefficientsQ12[coefficientIndex], LPC_QA - 12)
        coefficientIndex++
    }
    return lpcInversePredGainQA(workCoefficientsQA, order)
}

internal fun lpcInversePredGainQ24(
    predictionCoefficientsQ24: IntArray,
    order: Int
): Int {
    val workCoefficientsQA = arrayOf(IntArray(16), IntArray(16))
    val currentCoefficientsQA = workCoefficientsQA[order and 1]

    var coefficientIndex = 0
    while (coefficientIndex < order) {
        currentCoefficientsQA[coefficientIndex] =
            rshiftRound(predictionCoefficientsQ24[coefficientIndex], 24 - LPC_QA)
        coefficientIndex++
    }
    return lpcInversePredGainQA(workCoefficientsQA, order).invGainQ30
}

// Ported from tsilk.
// Source: schur.ts

internal fun schur(
    reflectionCoefficientsQ15: IntArray,
    correlations: IntArray,
    order: Int
): Int {
    val workCorrelations = IntArray((MAX_SCRATCH_LPC_ORDER + 1) * 2)
    val leadingZeros = clz32(correlations[0])

    var lag = 0
    while (lag <= order) {
        val normalizedCorrelation = when {
            leadingZeros < 2 -> rightShift(correlations[lag], 1)
            leadingZeros > 2 -> leftShift(correlations[lag], leadingZeros - 2)
            else -> correlations[lag]
        }
        workCorrelations[lag * 2] = normalizedCorrelation
        workCorrelations[lag * 2 + 1] = normalizedCorrelation
        lag++
    }

    var coefficientOrder = 0
    while (coefficientOrder < order) {
        val divisor = max32(rightShift(workCorrelations[1], 15), 1)
        var reflectionQ15 =
            -div3216(workCorrelations[(coefficientOrder + 1) * 2], divisor)
        reflectionQ15 = sat16(reflectionQ15)
        reflectionCoefficientsQ15[coefficientOrder] = reflectionQ15

        var lagOffset = 0
        while (lagOffset < order - coefficientOrder) {
            val firstCorrelation =
                workCorrelations[(lagOffset + coefficientOrder + 1) * 2]
            val secondCorrelation = workCorrelations[lagOffset * 2 + 1]
            workCorrelations[(lagOffset + coefficientOrder + 1) * 2] = smlawb(
                firstCorrelation,
                leftShift(secondCorrelation, 1),
                reflectionQ15
            )
            workCorrelations[lagOffset * 2 + 1] = smlawb(
                secondCorrelation,
                leftShift(firstCorrelation, 1),
                reflectionQ15
            )
            lagOffset++
        }
        coefficientOrder++
    }
    return workCorrelations[1]
}

// Ported from tsilk.
// Source: schur64.ts

internal fun schur64(
    reflectionCoefficientsQ16: IntArray,
    correlations: IntArray,
    order: Int
): Int {
    val workCorrelations = IntArray((MAX_SCRATCH_LPC_ORDER + 1) * 2)
    if (correlations[0] <= 0) {
        reflectionCoefficientsQ16.fill(0, 0, order)
        return 0
    }

    var lag = 0
    while (lag <= order) {
        workCorrelations[lag * 2] = correlations[lag]
        workCorrelations[lag * 2 + 1] = correlations[lag]
        lag++
    }

    var coefficientOrder = 0
    while (coefficientOrder < order) {
        val reflectionQ31 = div32VarQ(
            -workCorrelations[(coefficientOrder + 1) * 2],
            workCorrelations[1],
            31
        )
        reflectionCoefficientsQ16[coefficientOrder] =
            rshiftRound(reflectionQ31, 15)

        var lagOffset = 0
        while (lagOffset < order - coefficientOrder) {
            val firstCorrelation =
                workCorrelations[(lagOffset + coefficientOrder + 1) * 2]
            val secondCorrelation = workCorrelations[lagOffset * 2 + 1]
            workCorrelations[(lagOffset + coefficientOrder + 1) * 2] =
                firstCorrelation +
                        smmul(leftShift(secondCorrelation, 1), reflectionQ31)
            workCorrelations[lagOffset * 2 + 1] =
                secondCorrelation +
                        smmul(leftShift(firstCorrelation, 1), reflectionQ31)
            lagOffset++
        }
        coefficientOrder++
    }
    return workCorrelations[1]
}
