package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH

internal class VoiceActivityWorkspace(
    val bands: Array<IntArray> = Array(VAD_N_BANDS) { IntArray(MAX_FRAME_LENGTH / 2) },
    val bandEnergies: IntArray = IntArray(VAD_N_BANDS),
    val energyToNoiseRatiosQ8: IntArray = IntArray(VAD_N_BANDS),
)
