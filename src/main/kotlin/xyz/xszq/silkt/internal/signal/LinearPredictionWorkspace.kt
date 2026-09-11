package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.util.RefInt

internal class LinearPredictionWorkspace(
    val coefficientsQ16: IntArray,
    val tailCoefficientsQ16: IntArray,
    val candidateCoefficientsQ12: IntArray,
    val interpolatedNlsfQ15: IntArray,
    val filterState: IntArray,
    val residual: IntArray,
    val firstHalfEnergyOutput: RefInt,
    val secondHalfEnergyOutput: RefInt,
    val firstHalfShiftOutput: RefInt,
    val secondHalfShiftOutput: RefInt,
    val residualEnergyOutput: RefInt,
    val residualEnergyQOutput: RefInt,
    val tailEnergyOutput: RefInt,
    val tailEnergyQOutput: RefInt,
)
