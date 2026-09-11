package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.entropy.*
import xyz.xszq.silkt.internal.fixedpoint.div32
import xyz.xszq.silkt.internal.fixedpoint.fixConst
import xyz.xszq.silkt.internal.fixedpoint.limitInt
import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.quantization.processGains
import xyz.xszq.silkt.internal.quantization.quantizeNoiseShape
import xyz.xszq.silkt.internal.quantization.quantizeNoiseShapeDelayedDecision
import xyz.xszq.silkt.internal.signal.*
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.frameTerminationCdf
import xyz.xszq.silkt.internal.util.RefInt

internal const val SILK_LAST_FRAME: Int = 0
internal const val SILK_MORE_FRAMES: Int = 1

private const val SILK_LBRR__VER_1: Int = 2
private const val SILK_LBRR__VER_2: Int = 3
private const val SPEECH_ACTIVITY_DTX_THRESHOLD: Double = 0.1
private const val NO_SPEECH_FRAMES_BEFORE_DTX: Int = 5
private const val MAX_CONSECUTIVE_DTX: Int = 20
private const val SILK_ENC_PAYLOAD_BUF_TOO_SHORT: Int = -1
private const val SILK_ENC_INTERNAL_ERROR: Int = -2

internal fun encodeFrame(
    encoder: EncoderState,
    output: IntArray,
    encodedLength: RefInt,
    input: IntArray,
): Int {
    val runtime = encoder.runtime
    val analysis = encoder.analysis
    val frameCounterStart = runtime.frameCounter
    val inputState = encoder.input
    val workspace = inputState.frameWorkspace
    workspace.reset()
    val features = workspace.features
    val inputFrameOffset = encoder.frameGeometry.frameLength
    val residualPitch = workspace.residualPitch
    val residualPitchFrameOffset = encoder.frameGeometry.frameLength
    val prefilteredInput = workspace.prefilteredInput
    val highPassInput = workspace.highPassInput
    val speechActivityOutput = workspace.speechActivityOutput
    val snrDbOutput = workspace.snrDbOutput
    val tiltOutput = workspace.tiltOutput

    var status = getSpeechActivityQ8(
        analysis.voiceActivity,
        workspace.voiceActivity,
        speechActivityOutput,
        snrDbOutput,
        features.inputQualityBandsQ15,
        tiltOutput,
        input,
        encoder.frameGeometry.frameLength,
    )
    runtime.speechActivityQ8 = speechActivityOutput.value
    features.currentSnrDbQ7 = snrDbOutput.value
    features.inputTiltQ15 = tiltOutput.value
    features.seed = frameCounterStart and 3
    runtime.frameCounter = frameCounterStart + 1

    highPassVariableCutoff(
        analysis.noiseShaping,
        analysis.pitch,
        features,
        highPassInput,
        input,
        runtime.speechActivityQ8,
    )
    val laShapeOffset = inputFrameOffset + LA_SHAPE_MS * encoder.frameGeometry.samplingRateKHz
    lowPassVariableCutoff(
        inputState.transitionLowPassFilter,
        inputState.signalHistory,
        laShapeOffset,
        highPassInput,
        0,
        encoder.frameGeometry.frameLength,
    )
    findPitchLags(
        analysis.pitch,
        features,
        runtime.speechActivityQ8,
        residualPitch,
        inputState.signalHistory,
        inputFrameOffset,
    )
    noiseShapeAnalysis(
        analysis.noiseShaping,
        analysis.pitch,
        features,
        residualPitch,
        residualPitchFrameOffset,
        inputState.signalHistory,
        inputFrameOffset,
        runtime.speechActivityQ8,
        encoder.config.targetSnrDbQ7,
        runtime.channelBufferedMs,
        encoder.config.inBandFecSnrCompensationQ8,
    )
    prefilter(analysis.noiseShaping, features, prefilteredInput, inputState.signalHistory, inputFrameOffset)
    findPredictionCoefficients(
        analysis.prediction,
        features,
        residualPitch,
        inputState.signalHistory,
        runtime.firstFrameAfterReset,
        runtime.speechActivityQ8,
        encoder.config.packetLossPercentage,
        encoder.payload.frameCount,
        encoder.config.packetSizeMs,
    )
    processGains(
        analysis.noiseShaping.shape,
        encoder.frameGeometry.subframeLength,
        encoder.payload.frameCount,
        encoder.config.delayedDecisionCount,
        runtime.speechActivityQ8,
        features,
    )

    val lbrr = lbrrEncode(encoder, workspace, prefilteredInput)
    quantizeFrame(encoder, workspace, features, prefilteredInput)
    updateDtxState(encoder)

    if (encoder.payload.frameCount == 0) {
        rangeEncInit(encoder.payload.rangeCoder)
        encoder.payload.byteCount = 0
    }
    encodeParameters(
        encoder.payload.frameCount,
        encoder.frameGeometry.samplingRateKHz,
        runtime.previousTypeOffset,
        encoder.config.nlsfCodebooks,
        encoder.frameGeometry.frameLength,
        runtime.voiceActivityFlag,
        features,
        encoder.payload.rangeCoder,
        workspace.pulses,
        workspace.pulseEncoding,
    )
    advanceHistory(encoder, features)
    advancePayloadFrameCount(encoder)

    val payloadResult = finishPayload(encoder, workspace, output, encodedLength, features, lbrr)
    if (encoder.payload.rangeCoder.error != 0) {
        payloadResult.status = SILK_ENC_INTERNAL_ERROR
    }
    updateChannelBuffering(encoder, payloadResult.bytes)

    if (payloadResult.status != 0) {
        status = payloadResult.status
    }
    return status
}

private fun quantizeFrame(
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
            workspace.pulses,
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
        )
    } else {
        quantizeNoiseShape(
            noiseShapeQuantization,
            features,
            prefilteredInput,
            workspace.pulses,
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
        )
    }
}

private fun updateDtxState(encoder: EncoderState) {
    val runtime = encoder.runtime
    if (runtime.speechActivityQ8 < fixConst(SPEECH_ACTIVITY_DTX_THRESHOLD, 8)) {
        runtime.voiceActivityFlag = false
        runtime.silentFrameCount++
        if (runtime.silentFrameCount > NO_SPEECH_FRAMES_BEFORE_DTX) {
            runtime.inDtx = true
        }
        if (runtime.silentFrameCount > MAX_CONSECUTIVE_DTX + NO_SPEECH_FRAMES_BEFORE_DTX) {
            runtime.silentFrameCount = NO_SPEECH_FRAMES_BEFORE_DTX
            runtime.inDtx = false
        }
    } else {
        runtime.silentFrameCount = 0
        runtime.inDtx = false
        runtime.voiceActivityFlag = true
    }
}

private fun advanceHistory(encoder: EncoderState, features: FrameFeatures) {
    val moveSize = encoder.frameGeometry.frameLength + LA_SHAPE_MS * encoder.frameGeometry.samplingRateKHz
    encoder.input.signalHistory.copyInto(
        encoder.input.signalHistory,
        0,
        encoder.frameGeometry.frameLength,
        encoder.frameGeometry.frameLength + moveSize,
    )
    encoder.analysis.pitch.previousSignalType = features.signalType
    encoder.analysis.pitch.previousLag = features.pitchL[NB_SUBFR - 1]
    encoder.runtime.firstFrameAfterReset = false
}

private fun advancePayloadFrameCount(encoder: EncoderState) {
    if (encoder.payload.rangeCoder.error != 0) {
        encoder.payload.frameCount = 0
    } else {
        encoder.payload.frameCount++
    }
}

private fun finishPayload(
    encoder: EncoderState,
    workspace: FrameWorkspace,
    output: IntArray,
    encodedLength: RefInt,
    features: FrameFeatures,
    lbrr: EncodedPayload,
): PayloadResult {
    val packetDurationMs = encoder.payload.frameCount * FRAME_LENGTH_MS
    if (packetDurationMs < encoder.config.packetSizeMs) {
        return bufferPartialPayload(encoder, workspace, encodedLength)
    }

    val frameTerminator = selectFrameTerminator(encoder)
    rangeEncode(encoder.payload.rangeCoder, frameTerminator.terminator, frameTerminationCdf)
    val encodedBytes = workspace.encodedByteCount
    encodedBytes.value = 0
    rangeCoderGetLength(encoder.payload.rangeCoder, encodedBytes)
    var byteCount = encodedBytes.value
    val lbrrPackets = encoder.lbrr.packets

    if (encodedLength.value >= byteCount) {
        rangeEncWrapUp(encoder.payload.rangeCoder)
        for (index in 0 until byteCount) {
            output[index] = encoder.payload.rangeCoder.buffer[index]
        }
        if (frameTerminator.terminator > SILK_MORE_FRAMES &&
            lbrrPackets[frameTerminator.lbrrIndex].byteCount > 0 &&
            encodedLength.value >= byteCount + lbrrPackets[frameTerminator.lbrrIndex].byteCount
        ) {
            val extraBytes = lbrrPackets[frameTerminator.lbrrIndex].byteCount
            lbrrPackets[frameTerminator.lbrrIndex].payload.copyInto(output, byteCount, 0, extraBytes)
            byteCount += extraBytes
        }
        encodedLength.value = byteCount
        storeLbrrPacket(encoder, lbrrPackets, features, lbrr)
    } else {
        encodedLength.value = 0
        byteCount = 0
    }
    encoder.payload.frameCount = 0
    val status = if (byteCount == encodedLength.value) 0 else SILK_ENC_PAYLOAD_BUF_TOO_SHORT
    return PayloadResult(bytes = byteCount, status = status)
}

private fun selectFrameTerminator(encoder: EncoderState): FrameTerminator {
    val lbrrPackets = encoder.lbrr.packets
    var terminator = SILK_LAST_FRAME
    var lbrrIndex = (encoder.lbrr.oldestPacketIndex + 1) and LBRR_INDEX_MASK
    if (lbrrPackets[lbrrIndex].usage == SILK_ADD_LBRR_TO__PLUS_1) {
        terminator = SILK_LBRR__VER_1
    }
    if (lbrrPackets[encoder.lbrr.oldestPacketIndex].usage == SILK_ADD_LBRR_TO__PLUS_2) {
        terminator = SILK_LBRR__VER_2
        lbrrIndex = encoder.lbrr.oldestPacketIndex
    }
    return FrameTerminator(terminator = terminator, lbrrIndex = lbrrIndex)
}

private fun storeLbrrPacket(
    encoder: EncoderState,
    lbrrPackets: Array<LbrrPacket>,
    features: FrameFeatures,
    lbrr: EncodedPayload,
) {
    val index = encoder.lbrr.oldestPacketIndex
    lbrr.payload.copyInto(lbrrPackets[index].payload, 0, 0, lbrr.byteCount)
    lbrrPackets[index].byteCount = lbrr.byteCount
    lbrrPackets[index].usage = features.lbrrUsage
    encoder.lbrr.oldestPacketIndex = (encoder.lbrr.oldestPacketIndex + 1) and LBRR_INDEX_MASK
}

private fun bufferPartialPayload(
    encoder: EncoderState,
    workspace: FrameWorkspace,
    encodedLength: RefInt,
): PayloadResult {
    encodedLength.value = 0
    rangeEncode(encoder.payload.rangeCoder, SILK_MORE_FRAMES, frameTerminationCdf)
    val encodedBytes = workspace.encodedByteCount
    encodedBytes.value = 0
    rangeCoderGetLength(encoder.payload.rangeCoder, encodedBytes)
    return PayloadResult(bytes = encodedBytes.value, status = 0)
}

private fun updateChannelBuffering(encoder: EncoderState, byteCount: Int) {
    val runtime = encoder.runtime
    if (encoder.config.targetRateBps > 0) {
        runtime.channelBufferedMs += div32(
            8 * 1000 * (byteCount - encoder.payload.byteCount),
            encoder.config.targetRateBps,
        )
        runtime.channelBufferedMs -= FRAME_LENGTH_MS
        runtime.channelBufferedMs = limitInt(runtime.channelBufferedMs, 0, 100)
    }
    encoder.payload.byteCount = byteCount
}
