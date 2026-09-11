package xyz.xszq.silkt.internal.quantization

internal class DelayedDecisionState(
    var seedHistory: IntArray = IntArray(0),
    var pulseQ10: IntArray = IntArray(0),
    var quantizedSignalQ10: IntArray = IntArray(0),
    var lpcExcitationQ16: IntArray = IntArray(0),
    var ltpShapeQ10: IntArray = IntArray(0),
    var gainQ16: IntArray = IntArray(0),
    var autoregressiveStateQ14: IntArray = IntArray(0),
    var lpcStateQ14: IntArray = IntArray(0),
    var lowFrequencyShapingStateQ12: Int = 0,
    var seed: Int = 0,
    var frameStartSeed: Int = 0,
    var rateDistortionQ10: Int = 0,
)
