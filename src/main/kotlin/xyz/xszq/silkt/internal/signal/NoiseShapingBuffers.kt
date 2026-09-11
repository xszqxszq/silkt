package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.model.MAX_SHAPE_LPC_ORDER
import xyz.xszq.silkt.internal.model.MAX_SHAPE_LPC_WIN_LENGTH
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.util.RefInt

internal class NoiseShapingBuffers {
    val energyScale = RefInt()
    val residualEnergy = RefInt()
    val autocorrelation = IntArray(MAX_SHAPE_LPC_ORDER + 1)
    val reflectionCoefficientsQ16 = IntArray(MAX_SHAPE_LPC_ORDER)
    val ar1Q24 = IntArray(MAX_SHAPE_LPC_ORDER)
    val ar2Q24 = IntArray(MAX_SHAPE_LPC_ORDER)
    val windowedInput = IntArray(MAX_SHAPE_LPC_WIN_LENGTH)
    val warpedFilterStateQ14 = IntArray(MAX_SHAPE_LPC_ORDER + 1)
    val warpedCorrelationQ10 = LongArray(MAX_SHAPE_LPC_ORDER + 1)
    val harmonicFilteredQ12 = IntArray(MAX_FRAME_LENGTH / NB_SUBFR)
    val shapedResidual = IntArray(MAX_FRAME_LENGTH / NB_SUBFR + MAX_SHAPE_LPC_ORDER)
}
