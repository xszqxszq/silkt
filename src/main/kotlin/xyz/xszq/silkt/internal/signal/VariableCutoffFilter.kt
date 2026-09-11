package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.leftShift
import xyz.xszq.silkt.internal.fixedpoint.rightShift
import xyz.xszq.silkt.internal.fixedpoint.sat16
import xyz.xszq.silkt.internal.fixedpoint.smlawb
import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.transitionLpAQ28
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.transitionLpBQ28


// Ported from tsilk.
// Source: lp_variable_cutoff.ts

private const val NUMERATOR_TAPS = 3
private const val DENOMINATOR_TAPS = 2
private const val TRANSITION_SEGMENTS = 5
private const val TRANSITION_TIME_UP_MS = 5120
private const val TRANSITION_TIME_DOWN_MS = 2560
private const val TRANSITION_FRAMES_UP = TRANSITION_TIME_UP_MS / FRAME_LENGTH_MS
private const val TRANSITION_FRAMES_DOWN = TRANSITION_TIME_DOWN_MS / FRAME_LENGTH_MS

internal class TransitionLowPassState(
    val filterState: IntArray = IntArray(2),
    var transitionFrameNumber: Int = 0,
    var transitionMode: Int = 0,
)

@Suppress("DuplicatedCode")
private fun interpolateFilterTaps(
    numeratorQ28: IntArray,
    denominatorQ28: IntArray,
    segmentIndex: Int,
    interpolationQ16: Int,
) {
    if (segmentIndex >= TRANSITION_SEGMENTS - 1) {
        transitionLpBQ28[TRANSITION_SEGMENTS - 1].copyInto(numeratorQ28)
        transitionLpAQ28[TRANSITION_SEGMENTS - 1].copyInto(denominatorQ28)
        return
    }

    if (interpolationQ16 <= 0) {
        transitionLpBQ28[segmentIndex].copyInto(numeratorQ28)
        transitionLpAQ28[segmentIndex].copyInto(denominatorQ28)
        return
    }

    when (interpolationQ16) {
        sat16(interpolationQ16) -> {
            var tap = 0
            while (tap < NUMERATOR_TAPS) {
                val current = transitionLpBQ28[segmentIndex][tap]
                val next = transitionLpBQ28[segmentIndex + 1][tap]
                numeratorQ28[tap] = smlawb(current, next - current, interpolationQ16)
                tap++
            }
            tap = 0
            while (tap < DENOMINATOR_TAPS) {
                val current = transitionLpAQ28[segmentIndex][tap]
                val next = transitionLpAQ28[segmentIndex + 1][tap]
                denominatorQ28[tap] = smlawb(current, next - current, interpolationQ16)
                tap++
            }
        }

        1 shl 15 -> {
            var tap = 0
            while (tap < NUMERATOR_TAPS) {
                numeratorQ28[tap] = rightShift(
                    transitionLpBQ28[segmentIndex][tap] + transitionLpBQ28[segmentIndex + 1][tap],
                    1,
                )
                tap++
            }
            tap = 0
            while (tap < DENOMINATOR_TAPS) {
                denominatorQ28[tap] = rightShift(
                    transitionLpAQ28[segmentIndex][tap] + transitionLpAQ28[segmentIndex + 1][tap],
                    1,
                )
                tap++
            }
        }

        else -> {
            val inverseQ16 = (1 shl 16) - interpolationQ16
            var tap = 0
            while (tap < NUMERATOR_TAPS) {
                val current = transitionLpBQ28[segmentIndex][tap]
                val next = transitionLpBQ28[segmentIndex + 1][tap]
                numeratorQ28[tap] = smlawb(next, current - next, inverseQ16)
                tap++
            }
            tap = 0
            while (tap < DENOMINATOR_TAPS) {
                val current = transitionLpAQ28[segmentIndex][tap]
                val next = transitionLpAQ28[segmentIndex + 1][tap]
                denominatorQ28[tap] = smlawb(next, current - next, inverseQ16)
                tap++
            }
        }
    }
}

internal fun lowPassVariableCutoff(
    state: TransitionLowPassState,
    output: IntArray,
    outputOffset: Int,
    input: IntArray,
    inputOffset: Int,
    frameLength: Int,
) {
    val numeratorQ28 = IntArray(NUMERATOR_TAPS)
    val denominatorQ28 = IntArray(DENOMINATOR_TAPS)
    var interpolationQ16: Int
    var segmentIndex: Int

    if (state.transitionFrameNumber > 0) {
        if (state.transitionMode == 0) {
            if (state.transitionFrameNumber < TRANSITION_FRAMES_DOWN) {
                interpolationQ16 = leftShift(state.transitionFrameNumber, 16 - 5)
                segmentIndex = rightShift(interpolationQ16, 16)
                interpolationQ16 -= leftShift(segmentIndex, 16)
                interpolateFilterTaps(
                    numeratorQ28,
                    denominatorQ28,
                    segmentIndex,
                    interpolationQ16,
                )
                state.transitionFrameNumber++
            } else {
                interpolateFilterTaps(
                    numeratorQ28,
                    denominatorQ28,
                    TRANSITION_SEGMENTS - 1,
                    0,
                )
            }
        } else {
            if (state.transitionFrameNumber < TRANSITION_FRAMES_UP) {
                interpolationQ16 = leftShift(
                    TRANSITION_FRAMES_UP - state.transitionFrameNumber,
                    16 - 6,
                )
                segmentIndex = rightShift(interpolationQ16, 16)
                interpolationQ16 -= leftShift(segmentIndex, 16)
                interpolateFilterTaps(
                    numeratorQ28,
                    denominatorQ28,
                    segmentIndex,
                    interpolationQ16,
                )
                state.transitionFrameNumber++
            } else {
                interpolateFilterTaps(numeratorQ28, denominatorQ28, 0, 0)
            }
        }
    }

    if (state.transitionFrameNumber > 0) {
        applyBiquadFilter(
            input,
            inputOffset,
            numeratorQ28,
            denominatorQ28,
            state.filterState,
            0,
            output,
            outputOffset,
            frameLength,
        )
    } else {
        var sample = 0
        while (sample < frameLength) {
            output[outputOffset + sample] = input[inputOffset + sample]
            sample++
        }
    }
}
