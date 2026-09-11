package xyz.xszq.silkt.internal.quantization

internal class SampleState(
    @JvmField
    var pulseQ10: Int = 0,
    @JvmField
    var rateDistortionQ10: Int = 0,
    @JvmField
    var quantizedSignalQ14: Int = 0,
    @JvmField
    var lowFrequencyShapingStateQ12: Int = 0,
    @JvmField
    var ltpShapeQ10: Int = 0,
    @JvmField
    var lpcExcitationQ16: Int = 0,
)
