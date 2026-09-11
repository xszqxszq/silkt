package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.model.LTP_ORDER
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.model.MAX_SHAPE_LPC_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR

private const val NLSF_MSVQ_MAX_CB_STAGES: Int = 10

/** Per-frame analysis results shared by encoder stages. */
internal class FrameFeatures {
    val gainsQ16: IntArray = IntArray(NB_SUBFR)
    val predCoefQ12: IntArray = IntArray(2 * MAX_LPC_ORDER)
    val ltpCoefQ14: IntArray = IntArray(LTP_ORDER * NB_SUBFR)
    val residualEnergy: IntArray = IntArray(NB_SUBFR)
    val residualEnergyShift: IntArray = IntArray(NB_SUBFR)
    val ar1Q13: IntArray = IntArray(NB_SUBFR * MAX_SHAPE_LPC_ORDER)
    val ar2Q13: IntArray = IntArray(NB_SUBFR * MAX_SHAPE_LPC_ORDER)
    val lowFrequencyShapingQ14: IntArray = IntArray(NB_SUBFR)
    val gainsPreQ14: IntArray = IntArray(NB_SUBFR)
    val harmBoostQ14: IntArray = IntArray(NB_SUBFR)
    val tiltQ14: IntArray = IntArray(NB_SUBFR)
    val harmShapeGainQ14: IntArray = IntArray(NB_SUBFR)
    val inputQualityBandsQ15: IntArray = IntArray(4)
    val pitchL: IntArray = IntArray(NB_SUBFR)
    val ltpIndex: IntArray = IntArray(NB_SUBFR)
    val gainsIndices: IntArray = IntArray(NB_SUBFR)
    val nlsfIndices: IntArray = IntArray(NLSF_MSVQ_MAX_CB_STAGES)

    var ltpScaleQ14: Int = 0
    var ltPredCodGainQ7: Int = 0
    var nlsfInterpCoefQ2: Int = 0
    var lambdaQ10: Int = 0
    var inputQualityQ14: Int = 0
    var codingQualityQ14: Int = 0
    var pitchFreqLowHz: Int = 0
    var currentSnrDbQ7: Int = 0
    var signalType: Int = 0
    var quantOffsetType: Int = 0
    var lbrrUsage: Int = 0
    var seed: Int = 0
    var perIndex: Int = 0
    var ltpScaleIndex: Int = 0
    var lagIndex: Int = 0
    var contourIndex: Int = 0
    var predGainQ16: Int = 0
    var inputTiltQ15: Int = 0
    var sparsenessQ8: Int = 0

    @Suppress("DuplicatedCode")
    fun reset() {
        gainsQ16.fill(0)
        predCoefQ12.fill(0)
        ltpCoefQ14.fill(0)
        residualEnergy.fill(0)
        residualEnergyShift.fill(0)
        ar1Q13.fill(0)
        ar2Q13.fill(0)
        lowFrequencyShapingQ14.fill(0)
        gainsPreQ14.fill(0)
        harmBoostQ14.fill(0)
        tiltQ14.fill(0)
        harmShapeGainQ14.fill(0)
        inputQualityBandsQ15.fill(0)
        pitchL.fill(0)
        ltpIndex.fill(0)
        gainsIndices.fill(0)
        nlsfIndices.fill(0)

        ltpScaleQ14 = 0
        ltPredCodGainQ7 = 0
        nlsfInterpCoefQ2 = 0
        lambdaQ10 = 0
        inputQualityQ14 = 0
        codingQualityQ14 = 0
        pitchFreqLowHz = 0
        currentSnrDbQ7 = 0
        signalType = 0
        quantOffsetType = 0
        lbrrUsage = 0
        seed = 0
        perIndex = 0
        ltpScaleIndex = 0
        lagIndex = 0
        contourIndex = 0
        predGainQ16 = 0
        inputTiltQ15 = 0
        sparsenessQ8 = 0
    }
}
