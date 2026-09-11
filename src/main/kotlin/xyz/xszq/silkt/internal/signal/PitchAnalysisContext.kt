package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.model.SIG_TYPE_UNVOICED

internal class PitchAnalysisContext {
    val workspace = PitchAnalysisWorkspace()

    var frameLength: Int = 0
    var lookahead: Int = 0
    var sampleRateKHz: Int = 0
    var lpcOrder: Int = 16
    var complexity: Int = 2
    var correlationThresholdQ16: Int = 0
    var windowLength: Int = 0
    var previousLag: Int = 100
    var previousSignalType: Int = SIG_TYPE_UNVOICED
    var ltpCorrelationQ15: Int = 0

    fun configure(
        frameLength: Int,
        lookahead: Int,
        sampleRateKHz: Int,
        lpcOrder: Int,
        complexity: Int,
        correlationThresholdQ16: Int,
        windowLength: Int,
    ) {
        this.frameLength = frameLength
        this.lookahead = lookahead
        this.sampleRateKHz = sampleRateKHz
        this.lpcOrder = lpcOrder
        this.complexity = complexity
        this.correlationThresholdQ16 = correlationThresholdQ16
        this.windowLength = windowLength
    }

    fun reset() {
        previousLag = 100
        previousSignalType = SIG_TYPE_UNVOICED
        ltpCorrelationQ15 = 0
    }

    fun resetFrameHistory() {
        previousLag = 100
        previousSignalType = SIG_TYPE_UNVOICED
    }
}
