package xyz.xszq.silkt.internal.entropy

import xyz.xszq.silkt.internal.fixedpoint.clz32
import xyz.xszq.silkt.internal.fixedpoint.mulUint
import xyz.xszq.silkt.internal.model.MAX_ARITHM_BYTES
import xyz.xszq.silkt.internal.util.RefInt

private const val RANGE_CODER_WRITE_BEYOND_BUFFER: Int = -1
private const val RANGE_CODER_CDF_OUT_OF_RANGE: Int = -2
private const val RANGE_HIGH_BYTE_MASK: Int = -0x01000000
private const val RANGE_UPPER_TWO_BYTES_MASK: Int = -0x10000

/**
 * A negative index signals that carry has reached the byte before the output.
 * The encoder treats that boundary carry as fully propagated.
 */
private fun incrementByte(buffer: IntArray, index: Int): Int {
    if (index < 0) return 1
    val value = (buffer[index] + 1) and 0xFF
    buffer[index] = value
    return value
}

private fun isLessThanUnsigned(left: Int, right: Int): Boolean =
    (left xor Int.MIN_VALUE) < (right xor Int.MIN_VALUE)

// Ported from tsilk.
// Source: range_coder.ts

internal open class RangeCoderState {

    var bufferLength: Int = 0
    var bufferIndex: Int = 0
    var baseQ32: Int = 0
    var rangeQ16: Int = 0
    var error: Int = 0
    var buffer: IntArray = IntArray(MAX_ARITHM_BYTES)
}

internal fun rangeEncInit(rangeCoder: RangeCoderState) {
    rangeCoder.bufferLength = MAX_ARITHM_BYTES
    rangeCoder.rangeQ16 = 65535
    rangeCoder.bufferIndex = 0
    rangeCoder.baseQ32 = 0
    rangeCoder.error = 0
}

internal fun rangeEncode(
    rangeCoder: RangeCoderState,
    data: Int,
    prob: IntArray,
    probOffset: Int = 0
) {
    if (rangeCoder.error != 0) {
        return
    }
    if (probOffset < 0 || data < 0 || probOffset + data + 1 >= prob.size) {
        rangeCoder.error = RANGE_CODER_CDF_OUT_OF_RANGE
        return
    }
    val lowQ16 = prob[probOffset + data]
    val highQ16 = prob[probOffset + data + 1]
    val previousBaseQ32 = rangeCoder.baseQ32
    var baseQ32 = previousBaseQ32 + mulUint(rangeCoder.rangeQ16, lowQ16)
    val rangeQ32 = mulUint(rangeCoder.rangeQ16, highQ16 - lowQ16)
    if (isLessThanUnsigned(baseQ32, previousBaseQ32)) {
        var carryIndex = rangeCoder.bufferIndex
        while (true) {
            carryIndex--
            val incrementedByte = incrementByte(rangeCoder.buffer, carryIndex)
            if (incrementedByte != 0) {
                break
            }
        }
    }
    var rangeQ16: Int
    if (rangeQ32 and RANGE_HIGH_BYTE_MASK != 0) {
        rangeQ16 = rangeQ32 ushr 16
    } else {
        if (rangeQ32 and RANGE_UPPER_TWO_BYTES_MASK != 0) {
            rangeQ16 = rangeQ32 ushr 8
        } else {
            rangeQ16 = rangeQ32
            if (rangeCoder.bufferIndex >= rangeCoder.bufferLength) {
                rangeCoder.error = RANGE_CODER_WRITE_BEYOND_BUFFER
                return
            }
            rangeCoder.buffer[rangeCoder.bufferIndex++] = (baseQ32 ushr 24) and 255
            baseQ32 = baseQ32 shl 8
        }
        if (rangeCoder.bufferIndex >= rangeCoder.bufferLength) {
            rangeCoder.error = RANGE_CODER_WRITE_BEYOND_BUFFER
            return
        }
        rangeCoder.buffer[rangeCoder.bufferIndex++] = (baseQ32 ushr 24) and 255
        baseQ32 = baseQ32 shl 8
    }
    rangeCoder.baseQ32 = baseQ32
    rangeCoder.rangeQ16 = rangeQ16
}

internal fun rangeCoderGetLength(rangeCoder: RangeCoderState, byteCount: RefInt): Int {
    val bitCount = rangeCoder.bufferIndex * 8 + clz32(rangeCoder.rangeQ16 - 1) - 14
    byteCount.value = (bitCount + 7) shr 3
    return bitCount
}

internal fun rangeEncWrapUp(rangeCoder: RangeCoderState) {
    var baseQ24 = rangeCoder.baseQ32 ushr 8
    val byteCount = RefInt(value = 0)
    val bitsInStream = rangeCoderGetLength(rangeCoder, byteCount)
    val bitsToStore = bitsInStream - rangeCoder.bufferIndex * 8
    baseQ24 += 8388608 ushr (bitsToStore - 1)
    baseQ24 = baseQ24 and (-1 shl (24 - bitsToStore))
    if (baseQ24 and 16777216 != 0) {
        var carryIndex = rangeCoder.bufferIndex
        while (true) {
            if (incrementByte(rangeCoder.buffer, --carryIndex) != 0) {
                break
            }
        }
    }
    if (rangeCoder.bufferIndex < rangeCoder.bufferLength) {
        rangeCoder.buffer[rangeCoder.bufferIndex++] = (baseQ24 ushr 16) and 255
        if (bitsToStore > 8) {
            if (rangeCoder.bufferIndex < rangeCoder.bufferLength) {
                rangeCoder.buffer[rangeCoder.bufferIndex++] = (baseQ24 ushr 8) and 255
            }
        }
    }
    if (bitsInStream and 7 != 0) {
        val mask = 255 shr (bitsInStream and 7)
        val lastByteIndex = byteCount.value - 1
        if (lastByteIndex < rangeCoder.bufferLength) {
            rangeCoder.buffer[lastByteIndex] = rangeCoder.buffer[lastByteIndex] or mask
        }
    }
}
