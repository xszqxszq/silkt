package xyz.xszq.silkt.internal.signal

internal class PitchSelection(
    val lag: Int,
    val codebookIndex: Int,
    val ltpCorrelationQ15: Int,
)
