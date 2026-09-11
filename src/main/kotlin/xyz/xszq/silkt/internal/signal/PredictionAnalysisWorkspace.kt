package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.quantization.NlsfConversionWorkspace
import xyz.xszq.silkt.internal.quantization.NlsfQuantizationWorkspace

internal class PredictionAnalysisWorkspace(
    val inverseGainsQ16: IntArray,
    val localGains: IntArray,
    val subframeWeightsQ15: IntArray,
    val lpcInput: IntArray,
    val ltpWeights: IntArray,
    val ltpCoefficientsQ14: IntArray,
    val correlationShifts: IntArray,
    val nlsfQ15: IntArray,
    val nlsfWeightsQ6: IntArray,
    val interpolatedNlsfQ15: IntArray,
    val interpolatedWeightsQ6: IntArray,
    val nlsfIndices: IntArray,
    val firstPredictionCoefficientsQ12: IntArray,
    val secondPredictionCoefficientsQ12: IntArray,
    val predictionPairsQ12: Array<IntArray>,
    val predictionResidual: IntArray,
    val predictionFilterState: IntArray,
    val linearPrediction: LinearPredictionWorkspace,
    val longTermPrediction: LongTermPredictionWorkspace,
    val nlsfQuantization: NlsfQuantizationWorkspace,
    val nlsfConversion: NlsfConversionWorkspace,
)
