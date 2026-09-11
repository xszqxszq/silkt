package xyz.xszq.silkt.internal.signal

internal class NoiseShapingContext {
    var samplingRateKHz: Int = 16
    var frameLength: Int = 0
    var subframeLength: Int = 0
    var lookaheadShape: Int = 0
    var shapeWindowLength: Int = 0
    var shapingLpcOrder: Int = 16
    var warpingQ16: Int = 0
    val highPassFilterState: IntArray = IntArray(2)
    val highPassNumeratorQ28: IntArray = IntArray(3)
    val highPassDenominatorQ28: IntArray = IntArray(2)
    var highPassCutoffTrackingQ15: Int = 0
    var highPassCutoffSmoothedQ15: Int = 0
    var shape: ShapeState = ShapeState()
    var prefilter: PrefilterState = PrefilterState()
    var averageGainQ16: Int = 0
    val buffers = NoiseShapingBuffers()

    fun configure(
        samplingRateKHz: Int,
        frameLength: Int,
        subframeLength: Int,
        lookaheadShape: Int,
        shapeWindowLength: Int,
        shapingLpcOrder: Int,
        warpingQ16: Int,
    ) {
        this.samplingRateKHz = samplingRateKHz
        this.frameLength = frameLength
        this.subframeLength = subframeLength
        this.lookaheadShape = lookaheadShape
        this.shapeWindowLength = shapeWindowLength
        this.shapingLpcOrder = shapingLpcOrder
        this.warpingQ16 = warpingQ16
    }

    fun reset() {
        highPassFilterState.fill(0)
        highPassCutoffTrackingQ15 = 200844
        highPassCutoffSmoothedQ15 = 200844
        shape = ShapeState()
        prefilter = PrefilterState()
        averageGainQ16 = 0
    }

    fun resetFrameHistory() {
        prefilter.previousLag = 100
        shape.lastGainIndex = 1
    }
}
