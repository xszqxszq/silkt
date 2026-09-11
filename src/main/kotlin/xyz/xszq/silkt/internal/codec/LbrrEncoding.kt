package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.entropy.*
import xyz.xszq.silkt.internal.fixedpoint.fixConst
import xyz.xszq.silkt.internal.fixedpoint.limitInt
import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.model.MAX_ARITHM_BYTES
import xyz.xszq.silkt.internal.quantization.*
import xyz.xszq.silkt.internal.signal.FrameFeatures
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.frameTerminationCdf
import xyz.xszq.silkt.internal.util.RefInt

internal const val MAX_LBRR_DELAY: Int = 2
internal const val LBRR_LOSS_THRESHOLD: Int = 1
internal const val SILK_NO_LBRR: Int = 0
internal const val SILK_ADD_LBRR_TO__PLUS_1: Int = 1
internal const val SILK_ADD_LBRR_TO__PLUS_2: Int = 2

internal const val LBRR_INDEX_MASK: Int = MAX_LBRR_DELAY - 1

internal class LbrrPacket(
    var payload: IntArray,
    var byteCount: Int = 0,
    var usage: Int = 0,
)

private fun copyNoiseShapeQuantizerState(destination: NoiseShapeQuantizerState, source: NoiseShapeQuantizerState) {
    source.quantizedSignalQ10.copyInto(destination.quantizedSignalQ10)
    source.ltpShapeQ10.copyInto(destination.ltpShapeQ10)
    source.lpcStateQ14.copyInto(destination.lpcStateQ14)
    source.autoregressiveStateQ14.copyInto(destination.autoregressiveStateQ14)
    destination.lowFrequencyShapingStateQ12 = source.lowFrequencyShapingStateQ12
    destination.previousLag = source.previousLag
    destination.ltpExcitationWriteIndex = source.ltpExcitationWriteIndex
    destination.ltpShapeWriteIndex = source.ltpShapeWriteIndex
    destination.randomSeed = source.randomSeed
    destination.previousInverseGainQ16 = source.previousInverseGainQ16
    destination.shouldRewhiten = source.shouldRewhiten
}

private fun updateLbrrUsage(encoder: EncoderState, features: FrameFeatures) {
    if (encoder.config.lbrrEnabled) {
        var usage: Int = SILK_NO_LBRR
        if (encoder.runtime.speechActivityQ8 > fixConst(
                0.5,
                8
            ) && encoder.config.packetLossPercentage > LBRR_LOSS_THRESHOLD
        ) {
            usage = SILK_ADD_LBRR_TO__PLUS_1
        }
        features.lbrrUsage = usage
    } else {
        features.lbrrUsage = SILK_NO_LBRR
    }
}

private fun rateOnlyParameters(fsKHz: Int): Int = when (fsKHz) {
    8 -> 13500
    12 -> 15500
    16 -> 17500
    else -> 19500
}

internal fun lbrrEncode(
    encoder: EncoderState,
    workspace: FrameWorkspace,
    prefilteredInput: IntArray,
): EncodedPayload {
    val features = workspace.features
    val payload = workspace.lbrrPayload
    var encodedByteCount = 0
    updateLbrrUsage(encoder, features)
    if (!encoder.config.lbrrEnabled) {
        payload.byteCount = 0
        return payload
    }
    val tempGainsIndices = workspace.lbrrGainIndices
    val tempGainsQ16 = workspace.lbrrGainsQ16
    features.gainsIndices.copyInto(tempGainsIndices)
    features.gainsQ16.copyInto(tempGainsQ16)
    val typeOffsetPrevious = encoder.runtime.previousTypeOffset.value
    val ltpScaleIndex: Int = features.ltpScaleIndex
    if (
        encoder.config.complexity > 0 &&
        encoder.config.targetRateBps > rateOnlyParameters(encoder.frameGeometry.samplingRateKHz)
    ) {
        if (encoder.payload.frameCount == 0) {
            copyNoiseShapeQuantizerState(
                encoder.analysis.noiseShapeQuantization.lbrrState,
                encoder.analysis.noiseShapeQuantization.primaryState
            )
            encoder.lbrr.previousLastGainIndex = encoder.analysis.noiseShaping.shape.lastGainIndex
            features.gainsIndices[0] = limitInt(
                features.gainsIndices[0] + encoder.config.lbrrGainStepCount,
                0,
                GAIN_LEVEL_COUNT - 1,
            )
        }
        val previousLastGainIndex = RefInt(value = encoder.lbrr.previousLastGainIndex)
        dequantizeGains(
            features.gainsQ16,
            features.gainsIndices,
            previousLastGainIndex,
            if (encoder.payload.frameCount > 0) 1 else 0,
        )
        encoder.lbrr.previousLastGainIndex = previousLastGainIndex.value
        quantizeLbrrFrame(encoder, workspace, features, prefilteredInput)
    } else {
        workspace.lbrrPulses.fill(0)
        features.ltpScaleIndex = 0
    }
    if (encoder.payload.frameCount == 0) {
        rangeEncInit(encoder.lbrr.rangeCoder)
    }
    encodeParameters(
        encoder.payload.frameCount,
        encoder.frameGeometry.samplingRateKHz,
        encoder.runtime.previousTypeOffset,
        encoder.config.nlsfCodebooks,
        encoder.frameGeometry.frameLength,
        encoder.runtime.voiceActivityFlag,
        features,
        encoder.lbrr.rangeCoder,
        workspace.lbrrPulses,
        workspace.pulseEncoding,
    )
    val payloadFrameCount =
        if (encoder.lbrr.rangeCoder.error != 0) 0 else encoder.payload.frameCount + 1
    if (payloadFrameCount * FRAME_LENGTH_MS >= encoder.config.packetSizeMs) {
        rangeEncode(encoder.lbrr.rangeCoder, SILK_LAST_FRAME, frameTerminationCdf)
        val lenRet = workspace.encodedByteCount
        lenRet.value = 0
        rangeCoderGetLength(encoder.lbrr.rangeCoder, lenRet)
        if (lenRet.value <= MAX_ARITHM_BYTES) {
            rangeEncWrapUp(encoder.lbrr.rangeCoder)
            encoder.lbrr.rangeCoder.buffer.copyInto(payload.payload, 0, 0, lenRet.value)
            encodedByteCount = lenRet.value
        }
    } else {
        rangeEncode(encoder.lbrr.rangeCoder, SILK_MORE_FRAMES, frameTerminationCdf)
        encodedByteCount = 0
    }
    tempGainsIndices.copyInto(features.gainsIndices)
    tempGainsQ16.copyInto(features.gainsQ16)
    features.ltpScaleIndex = ltpScaleIndex
    encoder.runtime.previousTypeOffset.value = typeOffsetPrevious
    payload.byteCount = encodedByteCount
    return payload
}

private fun quantizeLbrrFrame(
    encoder: EncoderState,
    workspace: FrameWorkspace,
    features: FrameFeatures,
    prefilteredInput: IntArray,
) {
    val noiseShapeQuantization = encoder.analysis.noiseShapeQuantization
    if (noiseShapeQuantization.delayedDecisionCount > 1 || noiseShapeQuantization.warpingQ16 > 0) {
        quantizeNoiseShapeDelayedDecision(
            noiseShapeQuantization,
            features,
            prefilteredInput,
            workspace.lbrrPulses,
            features.nlsfInterpCoefQ2,
            features.predCoefQ12,
            features.ltpCoefQ14,
            features.ar2Q13,
            features.harmShapeGainQ14,
            features.tiltQ14,
            features.lowFrequencyShapingQ14,
            features.gainsQ16,
            features.lambdaQ10,
            features.ltpScaleQ14,
            true,
        )
    } else {
        quantizeNoiseShape(
            noiseShapeQuantization,
            features,
            prefilteredInput,
            workspace.lbrrPulses,
            features.nlsfInterpCoefQ2,
            features.predCoefQ12,
            features.ltpCoefQ14,
            features.ar2Q13,
            features.harmShapeGainQ14,
            features.tiltQ14,
            features.lowFrequencyShapingQ14,
            features.gainsQ16,
            features.lambdaQ10,
            features.ltpScaleQ14,
            true,
        )
    }
}
