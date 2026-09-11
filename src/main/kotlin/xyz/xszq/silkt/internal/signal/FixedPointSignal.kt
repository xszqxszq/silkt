package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*

// Ported from tsilk SILK sources.

internal fun movingAveragePrediction(
    input: IntArray,
    inputOffset: Int,
    bQ12: IntArray,
    bOffset: Int,
    state: IntArray,
    stateOffset: Int,
    output: IntArray,
    outputOffset: Int,
    length: Int,
    order: Int
) {
    var sampleIndex = 0
    while (sampleIndex < length) {
        val inputValue = input[inputOffset + sampleIndex]
        var outputValue = (inputValue shl 12) - state[stateOffset]
        outputValue = rshiftRound(outputValue, 12)
        var coefficientIndex = 0
        while (coefficientIndex < order - 1) {
            state[stateOffset + coefficientIndex] =
                smlabb(state[stateOffset + coefficientIndex + 1], inputValue, bQ12[bOffset + coefficientIndex])
            coefficientIndex++
        }
        state[stateOffset + order - 1] = smulbb(inputValue, bQ12[bOffset + order - 1])
        output[outputOffset + sampleIndex] = toInt16(sat16(outputValue))
        sampleIndex++
    }
}

internal fun expandBandwidth(
    coefficients: IntArray,
    coefficientOffset: Int,
    order: Int,
    chirpQ16: Int
) {
    val chirpMinusOneQ16 = chirpQ16 - 65536
    var chirp = chirpQ16
    var index = 0
    while (index < order - 1) {
        coefficients[coefficientOffset + index] =
            toInt16(rshiftRound(multiply(chirp, coefficients[coefficientOffset + index]), 16))
        chirp += rshiftRound(multiply(chirp, chirpMinusOneQ16), 16)
        index++
    }
    coefficients[coefficientOffset + order - 1] =
        toInt16(rshiftRound(multiply(chirp, coefficients[coefficientOffset + order - 1]), 16))
}

private val sigmoidSlopeQ10 = intArrayOf(237, 153, 73, 30, 12, 7)
private val sigmoidPositiveQ15 = intArrayOf(16384, 23955, 28861, 31213, 32178, 32548)
private val sigmoidNegativeQ15 = intArrayOf(16384, 8812, 3906, 1554, 589, 219)

internal fun sigmoidQ15(inputQ5: Int): Int {
    if (inputQ5 < 0) {
        val magnitude = -inputQ5
        if (magnitude >= 6 * 32) return 0
        val index = magnitude shr 5
        return sigmoidNegativeQ15[index] - smulbb(sigmoidSlopeQ10[index], magnitude and 31)
    }
    if (inputQ5 >= 6 * 32) return 32767
    val index = inputQ5 shr 5
    return sigmoidPositiveQ15[index] + smulbb(sigmoidSlopeQ10[index], inputQ5 and 31)
}

internal fun applyBiquadFilter(
    input: IntArray,
    inputOffset: Int,
    bQ28: IntArray,
    aQ28: IntArray,
    state: IntArray,
    stateOffset: Int,
    output: IntArray,
    outputOffset: Int,
    length: Int
) {
    val a0LowerQ28 = -aQ28[0] and 16383
    val a0UpperQ28 = -aQ28[0] shr 14
    val a1LowerQ28 = -aQ28[1] and 16383
    val a1UpperQ28 = -aQ28[1] shr 14
    var sampleIndex = 0
    while (sampleIndex < length) {
        val inputValue = input[inputOffset + sampleIndex]
        val outputQ14 = leftShift(smlawb(state[0], bQ28[0], inputValue), 2)
        state[stateOffset] =
            state[stateOffset + 1] + rshiftRound(smulwb(outputQ14, a0LowerQ28), 14)
        state[stateOffset] = smlawb(state[stateOffset], outputQ14, a0UpperQ28)
        state[stateOffset] = smlawb(state[stateOffset], bQ28[1], inputValue)
        state[stateOffset + 1] = rshiftRound(smulwb(outputQ14, a1LowerQ28), 14)
        state[stateOffset + 1] = smlawb(state[stateOffset + 1], outputQ14, a1UpperQ28)
        state[stateOffset + 1] = smlawb(state[stateOffset + 1], bQ28[2], inputValue)
        output[outputOffset + sampleIndex] = sat16((outputQ14 + (1 shl 14) - 1) shr 14)
        sampleIndex++
    }
}
