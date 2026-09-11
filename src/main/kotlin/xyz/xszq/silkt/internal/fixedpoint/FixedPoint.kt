package xyz.xszq.silkt.internal.fixedpoint

import xyz.xszq.silkt.internal.util.RefInt
import kotlin.math.pow

private const val int32Min: Int = Int.MIN_VALUE
private const val int16Max: Int = 32767
private const val int16Min: Int = -32768

internal const val INT32_MAX: Int = 2147483647

internal fun toInt32(value: Int): Int = value

internal fun toInt16(value: Int): Int = (value shl 16) shr 16

internal fun sat16(value: Int): Int = when {
    value > int16Max -> int16Max
    value < int16Min -> int16Min
    else -> value
}

internal fun rightShift(value: Int, shift: Int): Int = value shr shift

internal fun rshift32(value: Int, shift: Int): Int = rightShift(value, shift)

internal fun leftShift(value: Int, shift: Int): Int = toInt32(value shl shift)

internal fun rshiftRound(value: Int, shift: Int): Int {
    if (shift == 0) return value
    if (shift == 1) return (value shr 1) + (value and 1)
    return toInt32((value shr (shift - 1)) + 1) shr 1
}

internal fun multiply(a: Int, b: Int): Int = a * b

internal fun mulUint(a: Int, b: Int): Int {
    val product = a.toLong() * b.toLong()
    return (product and 0xffffffffL).toInt()
}

internal fun mla(a32: Int, b32: Int, c32: Int): Int = a32 + b32 * c32

internal fun smulwb(a32: Int, b16: Int): Int {
    val b32 = toInt16(b16)
    val highProduct = (a32 shr 16) * b32
    val lowProduct = toInt32((a32 and 65535) * b32) shr 16
    return highProduct + lowProduct
}

internal fun smlawb(a32: Int, b32: Int, c16: Int): Int =
    toInt32(a32 + smulwb(b32, c16))

internal fun smulwt(a32: Int, b32: Int): Int {
    val highProduct = (a32 shr 16) * (b32 shr 16)
    val lowProduct = toInt32((a32 and 65535) * (b32 shr 16)) shr 16
    return highProduct + lowProduct
}

internal fun smlawt(a32: Int, b32: Int, c32: Int): Int =
    toInt32(a32 + smulwt(b32, c32))

internal fun smulww(a32: Int, b32: Int): Int =
    mla(smulwb(a32, b32), a32, rshiftRound(b32, 16))

internal fun smmul(a32: Int, b32: Int): Int {
    val product = a32.toLong() * b32.toLong()
    return (product shr 32).toInt()
}

internal fun smulbb(a16: Int, b16: Int): Int =
    toInt32(toInt16(a16) * toInt16(b16))

internal fun div32(a: Int, b: Int): Int = a / b

internal fun div32VarQ(a: Int, b: Int, resultFractionBits: Int): Int {
    if (b == 0) return if (a >= 0) INT32_MAX else int32Min
    if (resultFractionBits < 0) return 0

    val absA = abs(a)
    val absB = abs(b)
    if (absA == 0) return 0

    val aHeadroom = clz32(absA) - 1
    val bHeadroom = clz32(absB) - 1
    var normalizedA = leftShift(a, aHeadroom)
    val normalizedB = leftShift(b, bHeadroom)
    val inverseB = div3216(INT32_MAX shr 2, rightShift(normalizedB, 16))
    var result = smulwb(normalizedA, inverseB)
    normalizedA -= smmul(normalizedB, result) shl 3
    result = smlawb(result, normalizedA, inverseB)

    val resultShift = 29 + aHeadroom - bHeadroom - resultFractionBits
    return when {
        resultShift <= 0 -> lshiftSat32(result, -resultShift)
        resultShift < 32 -> rightShift(result, resultShift)
        else -> 0
    }
}

internal fun rand(seed: Int): Int = toInt32(mla(907633515, seed, 196314165))

internal fun min(a: Int, b: Int): Int = if (a < b) a else b

internal fun max(a: Int, b: Int): Int = if (a > b) a else b

internal fun min32(a: Int, b: Int): Int = if (a < b) a else b

internal fun max32(a: Int, b: Int): Int = if (a > b) a else b

internal fun minInt(a: Int, b: Int): Int = if (a < b) a else b

internal fun maxInt(a: Int, b: Int): Int = if (a > b) a else b

internal fun max16(a: Int, b: Int): Int = if (a > b) a else b

internal fun limitInt(value: Int, lower: Int, upper: Int): Int = when {
    value < lower -> lower
    value > upper -> upper
    else -> value
}

internal fun limit32(value: Int, lower: Int, upper: Int): Int = when {
    value < lower -> lower
    value > upper -> upper
    else -> value
}

internal fun clz16(value: Int): Int = clz32(value and 65535) - 16

internal fun inverse32VarQ(b32: Int, resultFractionBits: Int): Int {
    if (b32 == 0) return INT32_MAX
    if (resultFractionBits <= 0) return 0

    val bHeadroom = clz32(abs(b32)) - 1
    val normalizedB = leftShift(b32, bHeadroom)
    val inverseB = div3216(INT32_MAX shr 2, rightShift(normalizedB, 16))
    var result = leftShift(inverseB, 16)
    val errorQ32 = -smulwb(normalizedB, inverseB) shl 3
    result = smlaww(result, errorQ32, inverseB)

    val resultShift = 61 - bHeadroom - resultFractionBits
    return when {
        resultShift <= 0 -> lshiftSat32(result, -resultShift)
        resultShift < 32 -> rightShift(result, resultShift)
        else -> 0
    }
}

internal fun lin2Log(linearValue: Int): Int {
    if (linearValue <= 0) return 0

    val leadingZeros = clz32(linearValue)
    val fractionQ7 = toInt32(
        rightShift(leftShift(linearValue, leadingZeros) and 2147483647, 24)
    )
    return toInt32(
        ((31 - leadingZeros) shl 7) +
                smlawb(fractionQ7, fractionQ7 * (128 - fractionQ7), 179)
    )
}

internal fun log2Lin(logarithmQ7: Int): Int {
    if (logarithmQ7 < 0) return 0
    if (logarithmQ7 >= 3968) return INT32_MAX

    val result = toInt32(1 shl (logarithmQ7 shr 7))
    val fractionQ7 = logarithmQ7 and 127
    val fractionAdjustment = smlawb(
        fractionQ7,
        fractionQ7 * (128 - fractionQ7),
        -174
    )
    return if (logarithmQ7 < 2048) {
        toInt32(result + rightShift(multiply(result, fractionAdjustment), 7))
    } else {
        toInt32(result + multiply(rightShift(result, 7), fractionAdjustment))
    }
}

internal fun abs(value: Int): Int = if (value < 0) -value else value

internal fun add32(a: Int, b: Int): Int = toInt32(a + b)

private fun add32Ovflw(a: Int, b: Int): Int =
    a + b

internal fun addRshift(a: Int, b: Int, shift: Int): Int =
    add32(a, rightShift(b, shift))

internal fun addLshift32(a: Int, b: Int, shift: Int): Int =
    add32(a, leftShift(b, shift))

internal fun lshift32(value: Int, shift: Int): Int = leftShift(value, shift)

internal fun subSat32(a: Int, b: Int): Int = a - b

internal fun sub32(a: Int, b: Int): Int = toInt32(a - b)

internal fun div3216(a32: Int, b16: Int): Int = a32 / b16

internal fun smlabb(a32: Int, b16: Int, c16: Int): Int =
    a32 + toInt16(b16) * toInt16(c16)

internal fun smlabbOvflw(a32: Int, b16: Int, c16: Int): Int =
    add32Ovflw(a32, smulbb(b16, c16))

internal fun addPosSat32(a: Int, b: Int): Int {
    val sum = toInt32(a + b)
    return if (sum < 0 && a > 0 && b > 0) INT32_MAX else sum
}

private fun ror32(value: Int, shift: Int): Int {
    val rotation = shift and 31
    if (rotation == 0) return value
    return (value ushr rotation) or (value shl (32 - rotation))
}

private class ClzFraction(var leadingZeros: Int = 0, var fractionQ7: Int = 0)

private fun clzFrac(value: Int): ClzFraction {
    val leadingZeros = clz32(value)
    val fractionQ7 = ror32(value, 24 - leadingZeros) and 127
    return ClzFraction(leadingZeros = leadingZeros, fractionQ7 = fractionQ7)
}

internal fun sqrtApprox(value: Int): Int {
    if (value <= 0) return 0

    val clzResult = clzFrac(value)
    val leadingZeros = clzResult.leadingZeros
    val fractionQ7 = clzResult.fractionQ7
    var root = if ((leadingZeros and 1) != 0) 32768 else 46214
    root = root shr rightShift(leadingZeros, 1)
    root = smlawb(root, root, smulbb(213, fractionQ7))
    return root
}

internal fun smlaww(a32: Int, b32: Int, c32: Int): Int =
    mla(smlawb(a32, b32, c32), b32, rshiftRound(c32, 16))

internal fun fixConst(value: Double, fractionBits: Int): Int =
    (value * 2.0.pow(fractionBits) + 0.5).toInt()

internal fun sumSquaresShift(
    energy: RefInt,
    scale: RefInt,
    samples: IntArray,
    samplesOffset: Int,
    length: Int
) {
    var shift = 0
    var energyValue = 0
    var sampleIndex = 0
    val finalSampleIndex = length - 1

    while (sampleIndex < finalSampleIndex) {
        energyValue += samples[samplesOffset + sampleIndex] * samples[samplesOffset + sampleIndex]
        energyValue += samples[samplesOffset + sampleIndex + 1] * samples[samplesOffset + sampleIndex + 1]
        sampleIndex += 2
        if (energyValue < 0) {
            energyValue = energyValue ushr 2
            shift = 2
            break
        }
    }

    while (sampleIndex < finalSampleIndex) {
        var pairEnergy = samples[samplesOffset + sampleIndex] * samples[samplesOffset + sampleIndex]
        pairEnergy += samples[samplesOffset + sampleIndex + 1] * samples[samplesOffset + sampleIndex + 1]
        energyValue += pairEnergy ushr shift
        if (energyValue < 0) {
            energyValue = energyValue ushr 2
            shift += 2
        }
        sampleIndex += 2
    }

    if (sampleIndex == finalSampleIndex) {
        val sampleEnergy = samples[samplesOffset + sampleIndex] * samples[samplesOffset + sampleIndex]
        energyValue += sampleEnergy ushr shift
    }

    if ((energyValue and -1073741824) != 0) {
        energyValue = energyValue ushr 2
        shift += 2
    }
    energy.value = energyValue
    scale.value = shift
}

internal fun lshiftSat32(value: Int, shift: Int): Int {
    if (shift <= 0) return toInt32(value)
    if (shift >= 31) {
        return when {
            value > 0 -> INT32_MAX
            value < 0 -> int32Min
            else -> 0
        }
    }

    val shifted = value.toLong() shl shift
    return when {
        shifted > INT32_MAX.toLong() -> INT32_MAX
        shifted < int32Min.toLong() -> int32Min
        else -> shifted.toInt()
    }
}

internal fun addSat16(a: Int, b: Int): Int {
    val sum = a + b
    return when {
        sum > int16Max -> int16Max
        sum < int16Min -> int16Min
        else -> sum
    }
}

internal fun clz32(value: Int): Int {
    if (value == 0) return 32

    var remaining = value
    var count = 0
    while (remaining >= 0) {
        count++
        remaining = remaining shl 1
    }
    return count
}
