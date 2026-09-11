package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.fixedpoint.max
import xyz.xszq.silkt.internal.fixedpoint.min
import xyz.xszq.silkt.internal.model.FRAME_LENGTH_MS
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.signal.initializeVoiceActivity
import xyz.xszq.silkt.internal.tables.NlsfTables.lpc10NlsfCodebooks
import xyz.xszq.silkt.internal.tables.NlsfTables.lpc16NlsfCodebooks
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.snrTableQ1
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.targetRateTableMb
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.targetRateTableNb
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.targetRateTableSwb
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.targetRateTableWb
import kotlin.math.roundToInt

// Ported from tsilk.
// Source: encoder.ts

private const val MIN_TARGET_RATE_BPS: Int = 5000
private const val MAX_TARGET_RATE_BPS: Int = 100000
private const val INBAND_FEC_MIN_RATE_BPS: Int = 18000
private const val LA_PITCH_MS: Int = 2
private const val FIND_PITCH_LPC_WIN_MS: Int = 20 + (LA_PITCH_MS shl 1)
private const val FIND_PITCH_CORRELATION_THRESHOLD_HC_MODE: Int = 45875
private const val FIND_PITCH_CORRELATION_THRESHOLD_MC_MODE: Int = 49152
private const val FIND_PITCH_CORRELATION_THRESHOLD_LC_MODE: Int = 52428
private const val TUNING_MU_LTP_QUANT_NB: Double = 0.03
private const val TUNING_MU_LTP_QUANT_MB: Double = 0.025
private const val TUNING_MU_LTP_QUANT_WB: Double = 0.02
private const val TUNING_MU_LTP_QUANT_SWB: Double = 0.016

internal const val LA_SHAPE_MS: Int = 5

private val supportedApiSampleRates: IntArray = intArrayOf(8000, 12000, 16000, 24000)

private val supportedPacketSizesMs: IntArray = intArrayOf(20, 40, 60, 80, 100)

private fun normalizeApiSampleRate(sampleRate: Int): Int {
    if (sampleRate !in supportedApiSampleRates) {
        throw IllegalStateException(
            Error(
                "Unsupported SILK sample rate: ${sampleRate}. Expected one of ${
                    supportedApiSampleRates.joinToString(
                        ", "
                    )
                } Hz."
            )
        )
    }
    return sampleRate
}

private fun normalizePacketSize(packetSizeMs: Int): Int {
    if (packetSizeMs !in supportedPacketSizesMs) {
        throw IllegalStateException(
            Error(
                "Unsupported SILK packet size: ${packetSizeMs}. Expected one of ${
                    supportedPacketSizesMs.joinToString(
                        ", "
                    )
                } ms."
            )
        )
    }
    return packetSizeMs
}

private fun toQ8(value: Double): Int {
    return (value * 256).roundToInt()
}

private fun computeTargetSnrQ7(sampleRateKHz: Int, targetRateBps: Int): Int {
    val rateTable = when (sampleRateKHz) {
        8 -> targetRateTableNb
        12 -> targetRateTableMb
        24 -> targetRateTableSwb
        else -> targetRateTableWb
    }
    val clampedRate = max(MIN_TARGET_RATE_BPS, min(MAX_TARGET_RATE_BPS, targetRateBps))
    var index = 1
    while (index < rateTable.size) {
        if (clampedRate <= rateTable[index]) {
            val rateLow = rateTable[index - 1]
            val rateHigh = rateTable[index]
            if (rateHigh == rateLow) {
                return snrTableQ1[index - 1] shl 6
            }
            val fractionQ6 = ((clampedRate - rateLow) shl 6) / (rateHigh - rateLow)
            return (snrTableQ1[index - 1] shl 6) +
                    fractionQ6 * (snrTableQ1[index] - snrTableQ1[index - 1])
        }
        index++
    }
    return snrTableQ1[snrTableQ1.size - 1] shl 6
}

private fun selectInitialInternalSampleRateKHz(
    apiSampleRate: Int,
    maxInternalSampleRate: Int,
    targetRateBps: Int
): Int {
    val rateThreshold = when {
        targetRateBps >= 25000 -> 24
        targetRateBps >= 14000 -> 16
        targetRateBps >= 10000 -> 12
        else -> 8
    }
    return min(
        rateThreshold,
        min(apiSampleRate / 1000, maxInternalSampleRate / 1000)
    )
}

internal open class SilkEncoder {

    lateinit var state: EncoderState

    constructor() {
        this.state = EncoderState()
        this.initialize()
    }

    fun initialize() {
        this.state = EncoderState()
        this.state.frameGeometry.configure(16, LA_PITCH_MS, LA_SHAPE_MS)
        initializeVoiceActivity(this.state.analysis.voiceActivity)
        this.state.runtime.firstFrameAfterReset = true
        this.state.analysis.pitch.reset()
        this.state.analysis.noiseShaping.reset()
        this.state.analysis.noiseShapeQuantization.primaryState.previousInverseGainQ16 = 65536
        this.state.analysis.noiseShapeQuantization.lbrrState.previousInverseGainQ16 = 65536
        this.state.analysis.noiseShapeQuantization.primaryState.previousLag = 100
        this.state.analysis.noiseShapeQuantization.lbrrState.previousLag = 100
        this.configureNlsfCodebooks()
        this.setupComplexity(2)
    }

    private fun setupComplexity(complexity: Int) {
        val level = max(0, min(2, complexity))
        val state = this.state
        val analysis = state.analysis
        var pitchComplexity: Int
        var pitchThresholdQ16: Int
        var pitchLpcOrder: Int
        var interpolateNlsfs: Boolean
        var lowComplexityLtpQuantization: Boolean
        var nlsfSurvivorCount: Int
        var shapingLpcOrder: Int
        var lookaheadShape: Int

        state.config.complexity = level
        when (level) {
            0 -> {
                pitchComplexity = 0
                pitchThresholdQ16 = FIND_PITCH_CORRELATION_THRESHOLD_LC_MODE
                pitchLpcOrder = 6
                shapingLpcOrder = 8
                lookaheadShape = 3 * state.frameGeometry.samplingRateKHz
                interpolateNlsfs = false
                lowComplexityLtpQuantization = true
                nlsfSurvivorCount = 2
            }

            1 -> {
                pitchComplexity = 1
                pitchThresholdQ16 = FIND_PITCH_CORRELATION_THRESHOLD_MC_MODE
                pitchLpcOrder = 12
                shapingLpcOrder = 12
                lookaheadShape = 5 * state.frameGeometry.samplingRateKHz
                interpolateNlsfs = false
                lowComplexityLtpQuantization = false
                nlsfSurvivorCount = 4
            }

            else -> {
                pitchComplexity = 2
                pitchThresholdQ16 = FIND_PITCH_CORRELATION_THRESHOLD_HC_MODE
                pitchLpcOrder = 16
                shapingLpcOrder = 16
                lookaheadShape = 5 * state.frameGeometry.samplingRateKHz
                interpolateNlsfs = true
                lowComplexityLtpQuantization = false
                nlsfSurvivorCount = 16
            }
        }
        state.config.delayedDecisionCount = when (level) {
            0 -> 1
            1 -> 2
            else -> 4
        }
        val warpingQ16 = if (level >= 1) state.frameGeometry.samplingRateKHz * 983 else 0
        state.frameGeometry.shapeLookahead = lookaheadShape
        pitchLpcOrder = min(
            pitchLpcOrder,
            if (state.config.predictionLpcOrder != 0) state.config.predictionLpcOrder else MAX_LPC_ORDER
        )
        val shapeWindowLength = 5 * state.frameGeometry.samplingRateKHz + 2 * lookaheadShape
        analysis.noiseShaping.configure(
            state.frameGeometry.samplingRateKHz,
            state.frameGeometry.frameLength,
            state.frameGeometry.subframeLength,
            lookaheadShape,
            shapeWindowLength,
            shapingLpcOrder,
            warpingQ16,
        )
        analysis.pitch.configure(
            state.frameGeometry.frameLength,
            state.frameGeometry.pitchLookahead,
            state.frameGeometry.samplingRateKHz,
            pitchLpcOrder,
            pitchComplexity,
            pitchThresholdQ16,
            FIND_PITCH_LPC_WIN_MS * state.frameGeometry.samplingRateKHz,
        )
        analysis.noiseShapeQuantization.configure(
            state.frameGeometry.frameLength,
            state.frameGeometry.subframeLength,
            state.config.predictionLpcOrder,
            shapingLpcOrder,
            state.config.delayedDecisionCount,
            warpingQ16,
        )
        analysis.prediction.configure(
            state.frameGeometry.frameLength,
            state.frameGeometry.subframeLength,
            state.config.predictionLpcOrder,
            analysis.prediction.ltpMuQ8,
            lowComplexityLtpQuantization,
            interpolateNlsfs,
            nlsfSurvivorCount,
            state.config.nlsfCodebooks,
        )
    }

    private fun configureNlsfCodebooks() {
        val config = this.state.config
        val useLpc10 = this.state.frameGeometry.samplingRateKHz <= 8
        val lpcOrder = if (useLpc10) 10 else 16
        config.predictionLpcOrder = lpcOrder
        config.nlsfCodebooks = if (useLpc10) lpc10NlsfCodebooks else lpc16NlsfCodebooks
    }

    fun configure(options: EncoderOptions) {
        val config = this.state.config
        val analysis = this.state.analysis
        val sampleRate =
            normalizeApiSampleRate(if (options.apiSampleRate != 0) options.apiSampleRate else 16000)
        val packetSizeMs = normalizePacketSize(options.packetSize ?: FRAME_LENGTH_MS)
        val maxInternalSampleRate = normalizeApiSampleRate(options.maxInternalSampleRate ?: sampleRate)
        val targetRate = options.bitRate ?: config.targetRateBps
        val fsKHz = selectInitialInternalSampleRateKHz(sampleRate, maxInternalSampleRate, targetRate)
        val fsChanged = this.state.frameGeometry.samplingRateKHz != fsKHz
        if (fsChanged) {
            this.state.payload.frameCount = 0
            this.state.payload.byteCount = 0
            this.state.lbrr.resetHistory()
            analysis.pitch.resetFrameHistory()
            this.state.runtime.firstFrameAfterReset = true
            analysis.noiseShaping.resetFrameHistory()
            analysis.noiseShapeQuantization.primaryState.previousLag = 100
            analysis.noiseShapeQuantization.primaryState.previousInverseGainQ16 = 65536
            analysis.noiseShapeQuantization.lbrrState.previousInverseGainQ16 = 65536
        }
        this.state.frameGeometry.configure(fsKHz, LA_PITCH_MS, LA_SHAPE_MS)
        config.targetRateBps = targetRate
        config.packetSizeMs = packetSizeMs
        config.packetLossPercentage = options.packetLossPercentage ?: 0
        config.useInBandFec = options.useInBandFec == true
        config.useDtx = options.useDtx == true
        config.targetSnrDbQ7 = computeTargetSnrQ7(fsKHz, targetRate)
        if (config.useInBandFec) {
            var lbrrRateThreshold = INBAND_FEC_MIN_RATE_BPS
            when (fsKHz) {
                8 -> {
                    lbrrRateThreshold -= 9000
                }

                12 -> {
                    lbrrRateThreshold -= 6000
                }

                16 -> {
                    lbrrRateThreshold -= 3000
                }
            }
            val packetLossPercentage = config.packetLossPercentage
            config.lbrrEnabled = true
            if (config.targetRateBps >= lbrrRateThreshold) {
                config.lbrrGainStepCount = max(8 - (packetLossPercentage shr 1), 0)
                if (packetLossPercentage > LBRR_LOSS_THRESHOLD) {
                    config.inBandFecSnrCompensationQ8 =
                        (6 shl 8) - (config.lbrrGainStepCount shl 7)
                } else {
                    config.inBandFecSnrCompensationQ8 = 0
                    config.lbrrEnabled = false
                }
            } else {
                config.inBandFecSnrCompensationQ8 = 0
                config.lbrrEnabled = false
            }
        } else {
            config.lbrrEnabled = false
            config.lbrrGainStepCount = 0
            config.inBandFecSnrCompensationQ8 = 0
        }
        analysis.prediction.ltpMuQ8 = when (fsKHz) {
            24 -> toQ8(TUNING_MU_LTP_QUANT_SWB)
            16 -> toQ8(TUNING_MU_LTP_QUANT_WB)
            12 -> toQ8(TUNING_MU_LTP_QUANT_MB)
            else -> toQ8(TUNING_MU_LTP_QUANT_NB)
        }
        this.configureNlsfCodebooks()
        this.setupComplexity(options.complexity ?: config.complexity)
    }

}
