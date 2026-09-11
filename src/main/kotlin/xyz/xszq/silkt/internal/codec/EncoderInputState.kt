package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.signal.TransitionLowPassState

internal class EncoderInputState {
    val transitionLowPassFilter: TransitionLowPassState = TransitionLowPassState()
    val signalHistory: IntArray = IntArray(2 * MAX_FRAME_LENGTH + 120)
    val frameWorkspace: FrameWorkspace = FrameWorkspace()
}
