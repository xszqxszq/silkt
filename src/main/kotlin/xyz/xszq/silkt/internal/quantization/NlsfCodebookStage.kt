package xyz.xszq.silkt.internal.quantization

internal class NlsfCodebookStage(
    val vectorCount: Int,
    val vectorsQ15: IntArray,
    val ratesQ5: IntArray
)
