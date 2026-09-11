package xyz.xszq.silkt.internal.signal

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class VariableCutoffFilterTest {
    @Test
    fun writesTransitionOutputToDestinationSlice() {
        val input = IntArray(20) { if (it % 2 == 0) 1000 else -500 }
        val output = IntArray(30) { 12345 }
        val state = TransitionLowPassState(transitionFrameNumber = 1)

        lowPassVariableCutoff(
            state = state,
            output = output,
            outputOffset = 10,
            input = input,
            inputOffset = 0,
            frameLength = input.size,
        )

        assertEquals(2, state.transitionFrameNumber)
        assertTrue(output.copyOfRange(10, 30).any { it != 12345 })
        assertTrue(output.copyOfRange(0, 10).all { it == 12345 })
    }
}
