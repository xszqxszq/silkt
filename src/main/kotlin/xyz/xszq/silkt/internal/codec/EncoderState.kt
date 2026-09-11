package xyz.xszq.silkt.internal.codec


// Ported from tsilk SILK sources.

// Ported from tsilk.
// Source: structs.ts

internal open class EncoderState {

    val frameGeometry: EncoderFrameGeometry = EncoderFrameGeometry()
    val config: EncoderConfig = EncoderConfig()
    val payload: EncoderPayloadState = EncoderPayloadState()
    val lbrr: EncoderLbrrState = EncoderLbrrState()
    val runtime: EncoderRuntimeState = EncoderRuntimeState()
    val input: EncoderInputState = EncoderInputState()
    val analysis: EncoderAnalysisState = EncoderAnalysisState()
}
