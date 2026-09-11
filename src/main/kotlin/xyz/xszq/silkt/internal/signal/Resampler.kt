package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import kotlin.math.ceil
import kotlin.math.floor
import kotlin.math.roundToInt

// Ported from tsilk.
// Source: resampler.ts

private val RESAMPLER_3_4_COEFFICIENTS = intArrayOf(
    -18249, -12532, -97, 284, -495, 309, 10268, 20317, -94, 156, -48, -720, 5984, 18278, -45, -4, 237, -847,
    2540, 14662,
)

private val RESAMPLER_2_3_COEFFICIENTS = intArrayOf(
    -11891, -12486, 20, 211, -657, 688, 8423, 15911, -44, 197, -152, -653, 3855, 13015
)

private val RESAMPLER_1_2_COEFFICIENTS =
    intArrayOf(2415, -13101, 158, -295, -400, 1265, 4832, 7968)

private val RESAMPLER_1_3_COEFFICIENTS =
    intArrayOf(16643, -14000, -331, 19, 581, 1421, 2290, 2845)

private val UPSAMPLE_2_HQ_FIRST_COEFFICIENTS = intArrayOf(4280, 33727 - 65536)
private val UPSAMPLE_2_HQ_SECOND_COEFFICIENTS = intArrayOf(16295, 54015 - 65536)
private val UPSAMPLE_2_HQ_NOTCH = intArrayOf(7864, -3604, 13107, 28508)

private const val DOWN_FIR_ORDER = 12
private const val LOW_QUALITY_MAX_BATCH_SIZE = 480
private const val DOWN_2_3_FIR_ORDER = 4
private const val DOWN_1_3_FIR_ORDER = 6

private fun upsample2HighQuality(input: IntArray): IntArray {
    val state = IntArray(6)
    val output = IntArray(input.size * 2)

    var sampleIndex = 0
    while (sampleIndex < input.size) {
        val inputSampleQ10 = input[sampleIndex] shl 10

        var difference = sub32(inputSampleQ10, state[0])
        var correction = smulwb(difference, UPSAMPLE_2_HQ_FIRST_COEFFICIENTS[0])
        var firstOutput = add32(state[0], correction)
        state[0] = add32(inputSampleQ10, correction)

        difference = sub32(firstOutput, state[1])
        correction = smlawb(
            difference,
            difference,
            UPSAMPLE_2_HQ_FIRST_COEFFICIENTS[1]
        )
        var secondOutput = add32(state[1], correction)
        state[1] = add32(firstOutput, correction)

        secondOutput = smlawb(
            secondOutput,
            state[5],
            UPSAMPLE_2_HQ_NOTCH[2]
        )
        secondOutput = smlawb(
            secondOutput,
            state[4],
            UPSAMPLE_2_HQ_NOTCH[1]
        )
        firstOutput = smlawb(secondOutput, state[4], UPSAMPLE_2_HQ_NOTCH[0])
        state[5] = sub32(secondOutput, state[5])
        output[sampleIndex * 2] = sat16(
            rshift32(smlawb(256, firstOutput, UPSAMPLE_2_HQ_NOTCH[3]), 9)
        )

        difference = sub32(inputSampleQ10, state[2])
        correction = smulwb(difference, UPSAMPLE_2_HQ_SECOND_COEFFICIENTS[0])
        firstOutput = add32(state[2], correction)
        state[2] = add32(inputSampleQ10, correction)

        difference = sub32(firstOutput, state[3])
        correction = smlawb(
            difference,
            difference,
            UPSAMPLE_2_HQ_SECOND_COEFFICIENTS[1]
        )
        secondOutput = add32(state[3], correction)
        state[3] = add32(firstOutput, correction)

        secondOutput = smlawb(
            secondOutput,
            state[4],
            UPSAMPLE_2_HQ_NOTCH[2]
        )
        secondOutput = smlawb(
            secondOutput,
            state[5],
            UPSAMPLE_2_HQ_NOTCH[1]
        )
        firstOutput = smlawb(secondOutput, state[5], UPSAMPLE_2_HQ_NOTCH[0])
        state[4] = sub32(secondOutput, state[4])
        output[sampleIndex * 2 + 1] = sat16(
            rshift32(smlawb(256, firstOutput, UPSAMPLE_2_HQ_NOTCH[3]), 9)
        )
        sampleIndex++
    }
    return output
}

private fun inverseRatioQ16(inputRate: Int, outputRate: Int): Int {
    return ceil(inputRate * 65536 / outputRate.toDouble()).toInt()
}

private fun upsampleByInterpolation(
    input: IntArray,
    ratio: Double,
    outputLength: Int
): IntArray {
    val output = IntArray(outputLength)
    var outputIndex = 0
    while (outputIndex < outputLength) {
        val sourcePosition = outputIndex / ratio
        val sourceIndex = floor(sourcePosition).toInt()
        val fraction = sourcePosition - sourceIndex
        val startSample = input[min(sourceIndex, input.size - 1)]
        val endSample = input[min(sourceIndex + 1, input.size - 1)]
        output[outputIndex] =
            (startSample + fraction * (endSample - startSample)).roundToInt()
        outputIndex++
    }
    return output
}

private fun downsampleByAverage(
    input: IntArray,
    ratio: Double,
    outputLength: Int
): IntArray {
    val output = IntArray(outputLength)
    var outputIndex = 0
    while (outputIndex < outputLength) {
        val sourceStart = outputIndex / ratio
        val sourceEnd = (outputIndex + 1) / ratio
        val firstSourceIndex = floor(sourceStart).toInt()
        val lastSourceIndex = min(ceil(sourceEnd).toInt(), input.size)
        var sum = 0
        var sampleCount = 0

        var sourceIndex = firstSourceIndex
        while (sourceIndex < lastSourceIndex) {
            sum += input[sourceIndex]
            sampleCount++
            sourceIndex++
        }
        output[outputIndex] = (sum / sampleCount).toDouble().roundToInt()
        outputIndex++
    }
    return output
}

private fun resampleDownFir(
    input: IntArray,
    inputRate: Int,
    outputRate: Int,
    coefficients: IntArray,
    firFractions: Int
): IntArray {
    if (input.isEmpty()) {
        return IntArray(0)
    }

    val batchSize = max(1, inputRate / 100)
    val iirState = IntArray(2)
    val firState = IntArray(DOWN_FIR_ORDER)
    val firCoefficients = coefficients.copyOfRange(2, coefficients.size)
    val indexIncrementQ16 = inverseRatioQ16(inputRate, outputRate)
    val outputSamples = mutableListOf<Int>()

    var inputOffset = 0
    while (inputOffset < input.size) {
        val sampleCount = min(batchSize, input.size - inputOffset)
        val filteredSamples = IntArray(sampleCount + DOWN_FIR_ORDER)
        firState.copyInto(filteredSamples)
        resamplerAR2(
            iirState,
            0,
            filteredSamples,
            DOWN_FIR_ORDER,
            input,
            inputOffset,
            coefficients,
            0,
            sampleCount
        )

        val maximumIndexQ16 = sampleCount shl 16
        var indexQ16 = 0
        while (indexQ16 < maximumIndexQ16) {
            val sampleIndex = indexQ16 shr 16
            var resultQ6: Int
            if (firFractions == 1) {
                resultQ6 = smulwb(
                    filteredSamples[sampleIndex] + filteredSamples[sampleIndex + 11],
                    firCoefficients[0]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 1] +
                            filteredSamples[sampleIndex + 10],
                    firCoefficients[1]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 2] +
                            filteredSamples[sampleIndex + 9],
                    firCoefficients[2]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 3] +
                            filteredSamples[sampleIndex + 8],
                    firCoefficients[3]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 4] +
                            filteredSamples[sampleIndex + 7],
                    firCoefficients[4]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 5] +
                            filteredSamples[sampleIndex + 6],
                    firCoefficients[5]
                )
            } else {
                val fraction = indexQ16 - (sampleIndex shl 16)
                val interpolationIndex = fraction * firFractions shr 16
                val lowOffset = interpolationIndex * 6
                val highOffset = (firFractions - 1 - interpolationIndex) * 6

                resultQ6 = smulwb(
                    filteredSamples[sampleIndex],
                    firCoefficients[lowOffset]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 1],
                    firCoefficients[lowOffset + 1]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 2],
                    firCoefficients[lowOffset + 2]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 3],
                    firCoefficients[lowOffset + 3]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 4],
                    firCoefficients[lowOffset + 4]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 5],
                    firCoefficients[lowOffset + 5]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 11],
                    firCoefficients[highOffset]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 10],
                    firCoefficients[highOffset + 1]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 9],
                    firCoefficients[highOffset + 2]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 8],
                    firCoefficients[highOffset + 3]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 7],
                    firCoefficients[highOffset + 4]
                )
                resultQ6 = smlawb(
                    resultQ6,
                    filteredSamples[sampleIndex + 6],
                    firCoefficients[highOffset + 5]
                )
            }
            outputSamples.add(sat16(rshiftRound(resultQ6, 6)))
            indexQ16 += indexIncrementQ16
        }

        filteredSamples.copyInto(
            firState,
            0,
            sampleCount,
            sampleCount + DOWN_FIR_ORDER
        )
        inputOffset += sampleCount
    }
    return outputSamples.toIntArray()
}

internal fun resample(
    input: IntArray,
    inputRate: Int,
    outputRate: Int
): IntArray {
    if (inputRate == outputRate) {
        return input.copyOf()
    }
    if ((inputRate == 24000 && outputRate == 12000) ||
        (inputRate == 16000 && outputRate == 8000)
    ) {
        return resampleDownFir(
            input,
            inputRate,
            outputRate,
            RESAMPLER_1_2_COEFFICIENTS,
            1
        )
    }
    if ((inputRate == 24000 && outputRate == 16000) ||
        (inputRate == 12000 && outputRate == 8000)
    ) {
        return resampleDownFir(
            input,
            inputRate,
            outputRate,
            RESAMPLER_2_3_COEFFICIENTS,
            2
        )
    }
    if (inputRate == 16000 && outputRate == 12000) {
        return resampleDownFir(
            input,
            inputRate,
            outputRate,
            RESAMPLER_3_4_COEFFICIENTS,
            3
        )
    }
    if (inputRate == 24000 && outputRate == 8000) {
        return resampleDownFir(
            input,
            inputRate,
            outputRate,
            RESAMPLER_1_3_COEFFICIENTS,
            1
        )
    }
    if (outputRate == inputRate * 2) {
        return upsample2HighQuality(input)
    }

    val ratio = outputRate.toDouble() / inputRate
    val outputLength = (input.size * ratio).roundToInt()
    if (ratio > 1) {
        return upsampleByInterpolation(input, ratio, outputLength)
    }
    return downsampleByAverage(input, ratio, outputLength)
}

// Ported from tsilk.
// Source: resampler_down.ts

private const val RESAMPLER_DOWN_2_FIRST_Q16 = 9872
private const val RESAMPLER_DOWN_2_SECOND_Q16 = 39809 - 65536

private val RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS =
    intArrayOf(-2797, -6507, 4697, 10739, 1567, 8276)

private val RESAMPLER_1_3_LOW_QUALITY_COEFFICIENTS =
    intArrayOf(16777, -9792, 890, 1614, 2148)

@Suppress("SameParameterValue")
private fun resamplerAR2(
    state: IntArray,
    stateOffset: Int,
    outputQ8: IntArray,
    outputOffset: Int,
    input: IntArray,
    inputOffset: Int,
    coefficientsQ14: IntArray,
    coefficientOffset: Int,
    length: Int
) {
    var sampleIndex = 0
    while (sampleIndex < length) {
        var outputQ32 = state[stateOffset] +
                (input[inputOffset + sampleIndex] shl 8)
        outputQ8[outputOffset + sampleIndex] = outputQ32
        outputQ32 = leftShift(outputQ32, 2)
        state[stateOffset] = smlawb(
            state[stateOffset + 1],
            outputQ32,
            coefficientsQ14[coefficientOffset]
        )
        state[stateOffset + 1] = smulwb(
            outputQ32,
            coefficientsQ14[coefficientOffset + 1]
        )
        sampleIndex++
    }
}

internal fun resamplerDown2(
    state: IntArray,
    stateOffset: Int,
    output: IntArray,
    outputOffset: Int,
    input: IntArray,
    inputOffset: Int,
    inputLength: Int
) {
    val outputPairCount = inputLength shr 1
    var inputSampleQ10: Int
    var outputSample: Int
    var difference: Int
    var correction: Int

    var pairIndex = 0
    while (pairIndex < outputPairCount) {
        inputSampleQ10 = input[inputOffset + pairIndex * 2] shl 10
        difference = sub32(inputSampleQ10, state[stateOffset])
        correction = smlawb(
            difference,
            difference,
            RESAMPLER_DOWN_2_SECOND_Q16
        )
        outputSample = add32(state[stateOffset], correction)
        state[stateOffset] = add32(inputSampleQ10, correction)

        inputSampleQ10 = input[inputOffset + pairIndex * 2 + 1] shl 10
        difference = sub32(inputSampleQ10, state[stateOffset + 1])
        correction = smulwb(difference, RESAMPLER_DOWN_2_FIRST_Q16)
        outputSample = add32(outputSample, state[stateOffset + 1])
        outputSample = add32(outputSample, correction)
        state[stateOffset + 1] = add32(inputSampleQ10, correction)
        output[outputOffset + pairIndex] = sat16(rshiftRound(outputSample, 11))
        pairIndex++
    }
}

internal fun resamplerDown23(
    state: IntArray,
    stateOffset: Int,
    output: IntArray,
    outputOffset: Int,
    input: IntArray,
    inputOffset: Int,
    inputLength: Int
) {
    var remainingInputLength = inputLength
    var inputSampleCount: Int
    val buffer = IntArray(LOW_QUALITY_MAX_BATCH_SIZE + DOWN_2_3_FIR_ORDER)
    state.copyInto(
        buffer,
        0,
        stateOffset,
        stateOffset + DOWN_2_3_FIR_ORDER
    )
    var inputReadIndex = 0
    var outputWriteIndex = 0

    while (true) {
        inputSampleCount = min(remainingInputLength, LOW_QUALITY_MAX_BATCH_SIZE)
        resamplerAR2(
            state,
            stateOffset + DOWN_2_3_FIR_ORDER,
            buffer,
            DOWN_2_3_FIR_ORDER,
            input,
            inputOffset + inputReadIndex,
            RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS,
            0,
            inputSampleCount
        )

        var bufferReadOffset = 0
        var remainingSamples = inputSampleCount
        while (remainingSamples > 2) {
            var resultQ6 = smulwb(
                buffer[bufferReadOffset],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[2]
            )
            resultQ6 = smlawb(
                resultQ6,
                buffer[bufferReadOffset + 1],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[3]
            )
            resultQ6 = smlawb(
                resultQ6,
                buffer[bufferReadOffset + 2],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[5]
            )
            resultQ6 = smlawb(
                resultQ6,
                buffer[bufferReadOffset + 3],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[4]
            )
            output[outputOffset + outputWriteIndex++] =
                sat16(rshiftRound(resultQ6, 6))

            resultQ6 = smulwb(
                buffer[bufferReadOffset + 1],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[4]
            )
            resultQ6 = smlawb(
                resultQ6,
                buffer[bufferReadOffset + 2],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[5]
            )
            resultQ6 = smlawb(
                resultQ6,
                buffer[bufferReadOffset + 3],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[3]
            )
            resultQ6 = smlawb(
                resultQ6,
                buffer[bufferReadOffset + 4],
                RESAMPLER_2_3_LOW_QUALITY_COEFFICIENTS[2]
            )
            output[outputOffset + outputWriteIndex++] =
                sat16(rshiftRound(resultQ6, 6))

            bufferReadOffset += 3
            remainingSamples -= 3
        }

        inputReadIndex += inputSampleCount
        remainingInputLength -= inputSampleCount
        if (remainingInputLength > 0) {
            buffer.copyInto(
                buffer,
                0,
                inputSampleCount,
                inputSampleCount + DOWN_2_3_FIR_ORDER
            )
        } else {
            break
        }
    }
    buffer.copyInto(
        state,
        stateOffset,
        inputSampleCount,
        inputSampleCount + DOWN_2_3_FIR_ORDER
    )
}

internal fun resamplerDown3(
    state: IntArray,
    stateOffset: Int,
    output: IntArray,
    outputOffset: Int,
    input: IntArray,
    inputOffset: Int,
    inputLength: Int
) {
    var remainingInputLength = inputLength
    var inputSampleCount: Int
    val buffer = IntArray(LOW_QUALITY_MAX_BATCH_SIZE + DOWN_1_3_FIR_ORDER)
    state.copyInto(
        buffer,
        0,
        stateOffset,
        stateOffset + DOWN_1_3_FIR_ORDER
    )
    var inputReadIndex = 0
    var outputWriteIndex = 0

    while (true) {
        inputSampleCount = min(remainingInputLength, LOW_QUALITY_MAX_BATCH_SIZE)
        resamplerAR2(
            state,
            stateOffset + DOWN_1_3_FIR_ORDER,
            buffer,
            DOWN_1_3_FIR_ORDER,
            input,
            inputOffset + inputReadIndex,
            RESAMPLER_1_3_LOW_QUALITY_COEFFICIENTS,
            0,
            inputSampleCount
        )

        var bufferReadOffset = 0
        var remainingSamples = inputSampleCount
        while (remainingSamples > 2) {
            var resultQ6 = smulwb(
                add32(
                    buffer[bufferReadOffset],
                    buffer[bufferReadOffset + 5]
                ),
                RESAMPLER_1_3_LOW_QUALITY_COEFFICIENTS[2]
            )
            resultQ6 = smlawb(
                resultQ6,
                add32(
                    buffer[bufferReadOffset + 1],
                    buffer[bufferReadOffset + 4]
                ),
                RESAMPLER_1_3_LOW_QUALITY_COEFFICIENTS[3]
            )
            resultQ6 = smlawb(
                resultQ6,
                add32(
                    buffer[bufferReadOffset + 2],
                    buffer[bufferReadOffset + 3]
                ),
                RESAMPLER_1_3_LOW_QUALITY_COEFFICIENTS[4]
            )
            output[outputOffset + outputWriteIndex++] =
                sat16(rshiftRound(resultQ6, 6))

            bufferReadOffset += 3
            remainingSamples -= 3
        }

        inputReadIndex += inputSampleCount
        remainingInputLength -= inputSampleCount
        if (remainingInputLength > 0) {
            buffer.copyInto(
                buffer,
                0,
                inputSampleCount,
                inputSampleCount + DOWN_1_3_FIR_ORDER
            )
        } else {
            break
        }
    }
    buffer.copyInto(
        state,
        stateOffset,
        inputSampleCount,
        inputSampleCount + DOWN_1_3_FIR_ORDER
    )
}
