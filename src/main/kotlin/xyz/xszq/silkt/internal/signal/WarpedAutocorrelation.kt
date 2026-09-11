package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.clz32
import xyz.xszq.silkt.internal.fixedpoint.smlawb
import xyz.xszq.silkt.internal.util.RefInt

// Ported from tsilk.
// Source: warped_autocorrelation.ts

private const val CORRELATION_FRACTION_BITS: Int = 10
private const val STATE_FRACTION_BITS: Int = 14
private const val CORRELATION_SHIFT: Int = 2 * STATE_FRACTION_BITS - CORRELATION_FRACTION_BITS

private fun countLeadingZeros64(value: Long): Int {
    if (value <= 0L) return 64
    val upper32 = (value shr 32 and 0xffffffffL).toInt()
    if (upper32 == 0) {
        val lower32 = (value and 0xffffffffL).toInt()
        return 32 + if (lower32 == 0) 32 else clz32(lower32)
    }
    return clz32(upper32)
}

internal fun warpedAutocorrelation(
    filterStateQ14: IntArray,
    accumulatedCorrelationQ10: LongArray,
    correlation: IntArray,
    scale: RefInt,
    input: IntArray,
    warpingQ16: Int,
    length: Int,
    order: Int,
) {
    val warpingCoefficientQ16 = warpingQ16 and 65535
    filterStateQ14.fill(0)
    accumulatedCorrelationQ10.fill(0)

    var sampleIndex = 0
    while (sampleIndex < length) {
        var temporaryQ14 = input[sampleIndex] shl STATE_FRACTION_BITS
        var coefficientIndex = 0
        while (coefficientIndex < order) {
            val firstTemporaryQ14 = smlawb(
                filterStateQ14[coefficientIndex],
                filterStateQ14[coefficientIndex + 1] - temporaryQ14,
                warpingCoefficientQ16,
            )
            filterStateQ14[coefficientIndex] = temporaryQ14
            accumulatedCorrelationQ10[coefficientIndex] += (
                    temporaryQ14.toLong() * filterStateQ14[0].toLong()
                    ) shr CORRELATION_SHIFT

            val secondTemporaryQ14 = smlawb(
                filterStateQ14[coefficientIndex + 1],
                filterStateQ14[coefficientIndex + 2] - firstTemporaryQ14,
                warpingCoefficientQ16,
            )
            filterStateQ14[coefficientIndex + 1] = firstTemporaryQ14
            accumulatedCorrelationQ10[coefficientIndex + 1] += (
                    firstTemporaryQ14.toLong() * filterStateQ14[0].toLong()
                    ) shr CORRELATION_SHIFT
            temporaryQ14 = secondTemporaryQ14
            coefficientIndex += 2
        }
        filterStateQ14[order] = temporaryQ14
        accumulatedCorrelationQ10[order] += (
                temporaryQ14.toLong() * filterStateQ14[0].toLong()
                ) shr CORRELATION_SHIFT
        sampleIndex++
    }

    val absoluteFirstCorrelation = if (accumulatedCorrelationQ10[0] >= 0L) {
        accumulatedCorrelationQ10[0]
    } else {
        -accumulatedCorrelationQ10[0]
    }
    var leftShift = countLeadingZeros64(absoluteFirstCorrelation) - 35
    val minimumLeftShift = -12 - CORRELATION_FRACTION_BITS
    val maximumLeftShift = 30 - CORRELATION_FRACTION_BITS
    if (leftShift < minimumLeftShift) {
        leftShift = minimumLeftShift
    } else if (leftShift > maximumLeftShift) {
        leftShift = maximumLeftShift
    }
    scale.value = -(CORRELATION_FRACTION_BITS + leftShift)

    for (index in 0..order) {
        correlation[index] = if (leftShift >= 0) {
            (accumulatedCorrelationQ10[index] shl leftShift).toInt()
        } else {
            (accumulatedCorrelationQ10[index] shr -leftShift).toInt()
        }
    }
}
