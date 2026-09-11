package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.model.LTP_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.util.RefInt

internal class LongTermPredictionWorkspace(
    val fittedCoefficientsQ16: IntArray,
    val coefficientDeltasQ14: IntArray,
    val coefficientSumsQ14: IntArray,
    val filteredResidualEnergies: IntArray,
    val diagonalWeights: IntArray,
    val crossCorrelation: IntArray,
    val signalEnergies: IntArray,
    val quantizationStageIndices: IntArray,
    val analysisCoefficientsQ14: IntArray,
    val energyOutput: RefInt,
    val shiftOutput: RefInt,
    val correlationShiftOutput: RefInt,
    val regularizedEnergy: RefInt,
    val quantizationBestIndex: RefInt,
    val quantizationRateDistortionQ14: RefInt,
    val cholesky: CholeskyWorkspace,
) {
    constructor() : this(
        fittedCoefficientsQ16 = IntArray(LTP_ORDER),
        coefficientDeltasQ14 = IntArray(LTP_ORDER),
        coefficientSumsQ14 = IntArray(NB_SUBFR),
        filteredResidualEnergies = IntArray(NB_SUBFR),
        diagonalWeights = IntArray(NB_SUBFR),
        crossCorrelation = IntArray(LTP_ORDER),
        signalEnergies = IntArray(NB_SUBFR),
        quantizationStageIndices = IntArray(NB_SUBFR),
        analysisCoefficientsQ14 = IntArray(LTP_ORDER),
        energyOutput = RefInt(),
        shiftOutput = RefInt(),
        correlationShiftOutput = RefInt(),
        regularizedEnergy = RefInt(),
        quantizationBestIndex = RefInt(),
        quantizationRateDistortionQ14 = RefInt(),
        cholesky = CholeskyWorkspace(),
    )
}
