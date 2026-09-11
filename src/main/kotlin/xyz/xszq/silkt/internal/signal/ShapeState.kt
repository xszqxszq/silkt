package xyz.xszq.silkt.internal.signal

internal class ShapeState(
    var lastGainIndex: Int = 1,
    var harmonicBoostSmoothingQ16: Int = 0,
    var harmonicShapeGainSmoothingQ16: Int = 0,
    var tiltSmoothingQ16: Int = 0,
)
