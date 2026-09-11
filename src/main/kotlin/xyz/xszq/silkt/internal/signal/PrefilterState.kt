package xyz.xszq.silkt.internal.signal

internal class PrefilterState(
    var ltpShape: IntArray = IntArray(512),
    var ltpShapeWriteIndex: Int = 0,
    var lowFrequencyShapingStateQ12: Int = 0,
    var lowFrequencyMovingAverageShapingQ12: Int = 0,
    var autoRegressiveShapingState: IntArray = IntArray(17),
    var harmonicHighPassState: Int = 0,
    var previousLag: Int = 100,
)
