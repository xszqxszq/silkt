package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*

// Ported from tsilk.
// Source: apply_sine_window.ts

private val frequencyTableQ16: IntArray = intArrayOf(
    12111, 9804, 8235, 7100, 6239, 5565, 5022, 4575, 4202, 3885, 3612, 3375, 3167, 2984, 2820, 2674, 2542, 2422, 2313,
    2214, 2123, 2038, 1961, 1889, 1822, 1760, 1702,
)

internal fun applySineWindow(
    windowedSignal: IntArray,
    windowedSignalOffset: Int,
    signal: IntArray,
    signalOffset: Int,
    windowType: Int,
    length: Int,
) {
    val frequencyQ16 = frequencyTableQ16[(length shr 2) - 4]
    val cosineQ16 = smulwb(frequencyQ16, -frequencyQ16)
    var sine0Q16: Int
    var sine1Q16: Int
    if (windowType == 1) {
        sine0Q16 = 0
        sine1Q16 = frequencyQ16 + (length shr 3)
    } else {
        sine0Q16 = 1 shl 16
        sine1Q16 = (1 shl 16) + (cosineQ16 shr 1) + (length shr 4)
    }

    var sampleIndex = 0
    while (sampleIndex < length) {
        val evenSamples = (signal[signalOffset + sampleIndex] and 65535) or
                (signal[signalOffset + sampleIndex + 1] shl 16)
        windowedSignal[windowedSignalOffset + sampleIndex] =
            smulwb(rightShift(sine0Q16 + sine1Q16, 1), evenSamples)
        windowedSignal[windowedSignalOffset + sampleIndex + 1] = smulwt(sine1Q16, evenSamples)
        sine0Q16 = (smulwb(sine1Q16, cosineQ16) + leftShift(sine1Q16, 1) - sine0Q16) + 1
        sine0Q16 = min(sine0Q16, 1 shl 16)

        val oddSamples = (signal[signalOffset + sampleIndex + 2] and 65535) or
                (signal[signalOffset + sampleIndex + 3] shl 16)
        windowedSignal[windowedSignalOffset + sampleIndex + 2] =
            smulwb(rightShift(sine0Q16 + sine1Q16, 1), oddSamples)
        windowedSignal[windowedSignalOffset + sampleIndex + 3] = smulwt(sine0Q16, oddSamples)
        sine1Q16 = (smulwb(sine0Q16, cosineQ16) + leftShift(sine0Q16, 1) - sine1Q16)
        sine1Q16 = min(sine1Q16, 1 shl 16)
        sampleIndex += 4
    }
}
