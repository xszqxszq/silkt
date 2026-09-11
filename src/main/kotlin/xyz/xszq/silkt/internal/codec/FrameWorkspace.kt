package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.entropy.PulseEncodingWorkspace
import xyz.xszq.silkt.internal.model.MAX_ARITHM_BYTES
import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.signal.FrameFeatures
import xyz.xszq.silkt.internal.signal.VoiceActivityWorkspace
import xyz.xszq.silkt.internal.util.RefInt

/** Reusable storage for data whose lifetime is one encoder frame. */
internal class FrameWorkspace {
    val features = FrameFeatures()
    val residualPitch = IntArray(2 * MAX_FRAME_LENGTH + 120)
    val prefilteredInput = IntArray(MAX_FRAME_LENGTH)
    val highPassInput = IntArray(MAX_FRAME_LENGTH)
    val pulses = IntArray(MAX_FRAME_LENGTH)
    val lbrrPulses = IntArray(MAX_FRAME_LENGTH)
    val pulseEncoding = PulseEncodingWorkspace()
    val voiceActivity = VoiceActivityWorkspace()
    val lbrrPayload = EncodedPayload(IntArray(MAX_ARITHM_BYTES))
    val lbrrGainIndices = IntArray(NB_SUBFR)
    val lbrrGainsQ16 = IntArray(NB_SUBFR)
    val speechActivityOutput = RefInt()
    val snrDbOutput = RefInt()
    val tiltOutput = RefInt()
    val encodedByteCount = RefInt()

    // Signal buffers are fully rewritten by the stages that own them.
    fun reset() {
        features.reset()
    }
}
