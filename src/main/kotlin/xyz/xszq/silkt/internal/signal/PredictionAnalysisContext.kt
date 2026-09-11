package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.model.LTP_ORDER
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.quantization.NlsfCodebook
import xyz.xszq.silkt.internal.quantization.NlsfConversionWorkspace
import xyz.xszq.silkt.internal.quantization.NlsfQuantizationWorkspace
import xyz.xszq.silkt.internal.util.RefInt

internal class PredictionAnalysisContext {
    var frameLength: Int = 0
    var subframeLength: Int = 0
    var predictionLpcOrder: Int = 16
    var ltpMuQ8: Int = 0
    var lowComplexityLtpQuantization: Boolean = false
    var interpolatedNlsfs: Boolean = false
    var nlsfSurvivorCount: Int = 0
    var nlsfCodebooks: Array<NlsfCodebook> = arrayOf()
    val state: PredictionState = PredictionState()
    val longTermPredictionGains: LongTermPredictionGainState = LongTermPredictionGainState()
    private var workspace: PredictionAnalysisWorkspace? = null

    internal val reusableWorkspace: PredictionAnalysisWorkspace
        get() = checkNotNull(workspace)

    fun configure(
        frameLength: Int,
        subframeLength: Int,
        predictionLpcOrder: Int,
        ltpMuQ8: Int,
        lowComplexityLtpQuantization: Boolean,
        interpolatedNlsfs: Boolean,
        nlsfSurvivorCount: Int,
        nlsfCodebooks: Array<NlsfCodebook>,
    ) {
        val effectivePredictionLpcOrder =
            if (predictionLpcOrder != 0) predictionLpcOrder else MAX_LPC_ORDER
        val maximumStageCount = nlsfCodebooks.maxOf { it.stageCount }
        val firstPredictionCoefficients = IntArray(effectivePredictionLpcOrder)
        val secondPredictionCoefficients = IntArray(effectivePredictionLpcOrder)
        workspace = PredictionAnalysisWorkspace(
            inverseGainsQ16 = IntArray(NB_SUBFR),
            localGains = IntArray(NB_SUBFR),
            subframeWeightsQ15 = IntArray(NB_SUBFR),
            lpcInput = IntArray(NB_SUBFR * (subframeLength + effectivePredictionLpcOrder)),
            ltpWeights = IntArray(NB_SUBFR * LTP_ORDER * LTP_ORDER),
            ltpCoefficientsQ14 = IntArray(NB_SUBFR * LTP_ORDER),
            correlationShifts = IntArray(NB_SUBFR),
            nlsfQ15 = IntArray(effectivePredictionLpcOrder),
            nlsfWeightsQ6 = IntArray(effectivePredictionLpcOrder),
            interpolatedNlsfQ15 = IntArray(effectivePredictionLpcOrder),
            interpolatedWeightsQ6 = IntArray(effectivePredictionLpcOrder),
            nlsfIndices = IntArray(maximumStageCount),
            firstPredictionCoefficientsQ12 = firstPredictionCoefficients,
            secondPredictionCoefficientsQ12 = secondPredictionCoefficients,
            predictionPairsQ12 = arrayOf(firstPredictionCoefficients, secondPredictionCoefficients),
            predictionResidual = IntArray(
                (subframeLength + effectivePredictionLpcOrder) * NB_SUBFR
            ),
            predictionFilterState = IntArray(MAX_LPC_ORDER),
            longTermPrediction = LongTermPredictionWorkspace(),
            nlsfQuantization = NlsfQuantizationWorkspace(),
            nlsfConversion = NlsfConversionWorkspace(),
            linearPrediction = LinearPredictionWorkspace(
                coefficientsQ16 = IntArray(effectivePredictionLpcOrder),
                tailCoefficientsQ16 = IntArray(effectivePredictionLpcOrder),
                candidateCoefficientsQ12 = IntArray(effectivePredictionLpcOrder),
                interpolatedNlsfQ15 = IntArray(effectivePredictionLpcOrder),
                filterState = IntArray(effectivePredictionLpcOrder),
                residual = IntArray(2 * (subframeLength + effectivePredictionLpcOrder)),
                firstHalfEnergyOutput = RefInt(),
                secondHalfEnergyOutput = RefInt(),
                firstHalfShiftOutput = RefInt(),
                secondHalfShiftOutput = RefInt(),
                residualEnergyOutput = RefInt(),
                residualEnergyQOutput = RefInt(),
                tailEnergyOutput = RefInt(),
                tailEnergyQOutput = RefInt(),
            ),
        )
        this.frameLength = frameLength
        this.subframeLength = subframeLength
        this.predictionLpcOrder = predictionLpcOrder
        this.ltpMuQ8 = ltpMuQ8
        this.lowComplexityLtpQuantization = lowComplexityLtpQuantization
        this.interpolatedNlsfs = interpolatedNlsfs
        this.nlsfSurvivorCount = nlsfSurvivorCount
        this.nlsfCodebooks = nlsfCodebooks
    }
}
