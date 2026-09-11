package xyz.xszq.silkt.internal.quantization

internal class SampleState(
    var pulseQ10: Int = 0,
    var rateDistortionQ10: Int = 0,
    var quantizedSignalQ14: Int = 0,
    var lowFrequencyShapingStateQ12: Int = 0,
    var ltpShapeQ10: Int = 0,
    var lpcExcitationQ16: Int = 0,
)
