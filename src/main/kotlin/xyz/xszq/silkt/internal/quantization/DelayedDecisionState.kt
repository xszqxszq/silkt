package xyz.xszq.silkt.internal.quantization

internal class DelayedDecisionState(
    @JvmField
    var seedHistory: IntArray = IntArray(0),
    @JvmField
    var pulseQ10: IntArray = IntArray(0),
    @JvmField
    var quantizedSignalQ10: IntArray = IntArray(0),
    @JvmField
    var lpcExcitationQ16: IntArray = IntArray(0),
    @JvmField
    var ltpShapeQ10: IntArray = IntArray(0),
    @JvmField
    var gainQ16: IntArray = IntArray(0),
    @JvmField
    var autoregressiveStateQ14: IntArray = IntArray(0),
    @JvmField
    var lpcStateQ14: IntArray = IntArray(0),
    @JvmField
    var lowFrequencyShapingStateQ12: Int = 0,
    @JvmField
    var seed: Int = 0,
    @JvmField
    var frameStartSeed: Int = 0,
    @JvmField
    var rateDistortionQ10: Int = 0,
)
