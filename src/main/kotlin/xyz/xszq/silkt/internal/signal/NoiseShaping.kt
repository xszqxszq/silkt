package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.*

private const val TUNING_FIND_PITCH_WHITE_NOISE_FRACTION = 0.001
private const val TUNING_LBRR_SPEECH_ACTIVITY_THRESHOLD = 0.5
private const val TUNING_BG_SNR_DECR_DB = 4.0
private const val TUNING_HARM_SNR_INCR_DB = 2.0
private const val TUNING_SPARSE_SNR_INCR_DB = 2.0
private const val TUNING_SPARSENESS_THRESHOLD_QNT_OFFSET = 0.75
private const val TUNING_SHAPE_WHITE_NOISE_FRACTION = 0.00001
private const val TUNING_BANDWIDTH_EXPANSION = 0.95
private const val TUNING_LOW_RATE_BANDWIDTH_EXPANSION_DELTA = 0.01
private const val TUNING_DE_ESSER_COEF_SWB_DB = 2.0
private const val TUNING_DE_ESSER_COEF_WB_DB = 1.0
private const val TUNING_LOW_RATE_HARMONIC_BOOST = 0.1
private const val TUNING_LOW_INPUT_QUALITY_HARMONIC_BOOST = 0.1
private const val TUNING_HARMONIC_SHAPING = 0.3
private const val TUNING_HIGH_RATE_OR_LOW_QUALITY_HARMONIC_SHAPING = 0.2
private const val TUNING_HP_NOISE_COEF = 0.3
private const val TUNING_HARM_HP_NOISE_COEF = 0.35
private const val TUNING_INPUT_TILT = 0.05
private const val TUNING_HIGH_RATE_INPUT_TILT = 0.1
private const val TUNING_LOW_FREQ_SHAPING = 3.0
private const val TUNING_LOW_QUALITY_LOW_FREQ_SHAPING_DECR = 0.5
private const val TUNING_NOISE_FLOOR_DB = 4.0
private const val TUNING_RELATIVE_MIN_GAIN_DB = -50.0
private const val TUNING_GAIN_SMOOTHING_COEF = 0.001
private const val TUNING_SUBFR_SMTH_COEF = 0.4

// Ported from tsilk.
// Source: noise_shape_analysis_FIX.ts

private fun warpedGain(coefficientsQ24: IntArray, lambdaQ16: Int, order: Int): Int {
    val warpedLambdaQ16 = -lambdaQ16
    var gainQ24 = coefficientsQ24[order - 1]
    var coefficientIndex = order - 2
    while (coefficientIndex >= 0) {
        gainQ24 = smlawb(coefficientsQ24[coefficientIndex], gainQ24, warpedLambdaQ16)
        coefficientIndex--
    }
    gainQ24 = smlawb(fixConst(1.0, 24), gainQ24, lambdaQ16)
    return inverse32VarQ(gainQ24, 40)
}

@Suppress("DuplicatedCode")
private fun limitWarpedCoefficients(
    synthesisCoefficientsQ24: IntArray,
    analysisCoefficientsQ24: IntArray,
    lambdaQ16: Int,
    limitQ24: Int,
    order: Int
) {
    var warpedLambdaQ16 = -lambdaQ16
    var numeratorQ16: Int
    var synthesisGainQ16: Int
    var analysisGainQ16: Int

    warpCoefficientPairsBackward(
        synthesisCoefficientsQ24,
        analysisCoefficientsQ24,
        warpedLambdaQ16,
        order
    )

    warpedLambdaQ16 = -warpedLambdaQ16
    numeratorQ16 = smlawb(fixConst(1.0, 16), -warpedLambdaQ16, warpedLambdaQ16)
    synthesisGainQ16 = pairedWarpedGainQ16(
        synthesisCoefficientsQ24,
        warpedLambdaQ16,
        numeratorQ16
    )
    analysisGainQ16 = pairedWarpedGainQ16(
        analysisCoefficientsQ24,
        warpedLambdaQ16,
        numeratorQ16
    )
    scaleCoefficientPairs(
        synthesisCoefficientsQ24,
        analysisCoefficientsQ24,
        synthesisGainQ16,
        analysisGainQ16,
        order
    )

    var iterationIndex = 0
    while (iterationIndex < 10) {
        val maximumCoefficientIndex = findMaximumPairedMagnitudeIndex(
            synthesisCoefficientsQ24,
            analysisCoefficientsQ24,
            order
        )
        val maximumMagnitudeQ24 = max(
            abs(synthesisCoefficientsQ24[maximumCoefficientIndex]),
            abs(analysisCoefficientsQ24[maximumCoefficientIndex])
        )
        if (maximumMagnitudeQ24 <= limitQ24) {
            return
        }

        warpCoefficientPairsForward(
            synthesisCoefficientsQ24,
            analysisCoefficientsQ24,
            warpedLambdaQ16,
            order
        )

        synthesisGainQ16 = inverse32VarQ(synthesisGainQ16, 32)
        analysisGainQ16 = inverse32VarQ(analysisGainQ16, 32)
        scaleCoefficientPairs(
            synthesisCoefficientsQ24,
            analysisCoefficientsQ24,
            synthesisGainQ16,
            analysisGainQ16,
            order
        )

        val bandwidthChirpQ16 = fixConst(0.99, 16) - div32VarQ(
            smulwb(
                maximumMagnitudeQ24 - limitQ24,
                smlabb(fixConst(0.8, 10), fixConst(0.1, 10), iterationIndex)
            ),
            multiply(maximumMagnitudeQ24, maximumCoefficientIndex + 1),
            22
        )
        expandBandwidth32(synthesisCoefficientsQ24, order, bandwidthChirpQ16)
        expandBandwidth32(analysisCoefficientsQ24, order, bandwidthChirpQ16)

        warpedLambdaQ16 = -warpedLambdaQ16
        warpCoefficientPairsBackward(
            synthesisCoefficientsQ24,
            analysisCoefficientsQ24,
            warpedLambdaQ16,
            order
        )
        warpedLambdaQ16 = -warpedLambdaQ16

        numeratorQ16 = smlawb(fixConst(1.0, 16), -warpedLambdaQ16, warpedLambdaQ16)
        synthesisGainQ16 = pairedWarpedGainQ16(
            synthesisCoefficientsQ24,
            warpedLambdaQ16,
            numeratorQ16
        )
        analysisGainQ16 = pairedWarpedGainQ16(
            analysisCoefficientsQ24,
            warpedLambdaQ16,
            numeratorQ16
        )
        scaleCoefficientPairs(
            synthesisCoefficientsQ24,
            analysisCoefficientsQ24,
            synthesisGainQ16,
            analysisGainQ16,
            order
        )
        iterationIndex++
    }
}

private fun warpCoefficientPairsBackward(
    synthesisCoefficientsQ24: IntArray,
    analysisCoefficientsQ24: IntArray,
    warpedLambdaQ16: Int,
    order: Int
) {
    var coefficientIndex = order - 1
    while (coefficientIndex > 0) {
        synthesisCoefficientsQ24[coefficientIndex - 1] = smlawb(
            synthesisCoefficientsQ24[coefficientIndex - 1],
            synthesisCoefficientsQ24[coefficientIndex],
            warpedLambdaQ16
        )
        analysisCoefficientsQ24[coefficientIndex - 1] = smlawb(
            analysisCoefficientsQ24[coefficientIndex - 1],
            analysisCoefficientsQ24[coefficientIndex],
            warpedLambdaQ16
        )
        coefficientIndex--
    }
}

private fun warpCoefficientPairsForward(
    synthesisCoefficientsQ24: IntArray,
    analysisCoefficientsQ24: IntArray,
    warpedLambdaQ16: Int,
    order: Int
) {
    var coefficientIndex = 1
    while (coefficientIndex < order) {
        synthesisCoefficientsQ24[coefficientIndex - 1] = smlawb(
            synthesisCoefficientsQ24[coefficientIndex - 1],
            synthesisCoefficientsQ24[coefficientIndex],
            warpedLambdaQ16
        )
        analysisCoefficientsQ24[coefficientIndex - 1] = smlawb(
            analysisCoefficientsQ24[coefficientIndex - 1],
            analysisCoefficientsQ24[coefficientIndex],
            warpedLambdaQ16
        )
        coefficientIndex++
    }
}

private fun pairedWarpedGainQ16(
    coefficientsQ24: IntArray,
    warpedLambdaQ16: Int,
    numeratorQ16: Int
): Int {
    val denominatorQ24 = smlawb(fixConst(1.0, 24), coefficientsQ24[0], warpedLambdaQ16)
    return div32VarQ(numeratorQ16, denominatorQ24, 24)
}

private fun scaleCoefficientPairs(
    synthesisCoefficientsQ24: IntArray,
    analysisCoefficientsQ24: IntArray,
    synthesisGainQ16: Int,
    analysisGainQ16: Int,
    order: Int
) {
    var coefficientIndex = 0
    while (coefficientIndex < order) {
        synthesisCoefficientsQ24[coefficientIndex] =
            smulww(synthesisGainQ16, synthesisCoefficientsQ24[coefficientIndex])
        analysisCoefficientsQ24[coefficientIndex] =
            smulww(analysisGainQ16, analysisCoefficientsQ24[coefficientIndex])
        coefficientIndex++
    }
}

private fun findMaximumPairedMagnitudeIndex(
    synthesisCoefficientsQ24: IntArray,
    analysisCoefficientsQ24: IntArray,
    order: Int
): Int {
    var maximumCoefficientIndex = 0
    var maximumMagnitudeQ24 = -1
    var coefficientIndex = 0
    while (coefficientIndex < order) {
        val magnitudeQ24 = max(
            abs(synthesisCoefficientsQ24[coefficientIndex]),
            abs(analysisCoefficientsQ24[coefficientIndex])
        )
        if (magnitudeQ24 > maximumMagnitudeQ24) {
            maximumMagnitudeQ24 = magnitudeQ24
            maximumCoefficientIndex = coefficientIndex
        }
        coefficientIndex++
    }
    return maximumCoefficientIndex
}

internal fun noiseShapeAnalysis(
    context: NoiseShapingContext,
    pitchAnalysis: PitchAnalysisContext,
    features: FrameFeatures,
    pitchResidual: IntArray,
    pitchResidualOffset: Int,
    input: IntArray,
    inputOffset: Int,
    speechActivityQ8: Int,
    targetSnrQ7: Int,
    bufferedChannelMilliseconds: Int,
    inBandFecSnrCompQ8: Int
) {
    val shapeState = context.shape
    val snrAdjustmentDBQ7 = analyzeNoiseShapeQuality(
        context,
        pitchAnalysis,
        features,
        pitchResidual,
        pitchResidualOffset,
        speechActivityQ8,
        targetSnrQ7,
        bufferedChannelMilliseconds,
        inBandFecSnrCompQ8
    )
    analyzeShapingCoefficients(context, features, input, inputOffset)
    adjustQuantizationGains(context, features, speechActivityQ8, snrAdjustmentDBQ7)
    applyPreEmphasisGain(context, features, speechActivityQ8)
    val tiltQ16 = calculateShapingTilt(context, features, speechActivityQ8)
    smoothShapingControls(shapeState, features, pitchAnalysis, tiltQ16)
}

private fun adjustQuantizationGains(
    context: NoiseShapingContext,
    features: FrameFeatures,
    speechActivityQ8: Int,
    snrAdjustmentDBQ7: Int
) {
    val gainMultiplierQ16 = log2Lin(
        -smlawb(-fixConst(16.0, 7), snrAdjustmentDBQ7, fixConst(0.16, 16))
    )
    var gainOffsetQ16 = log2Lin(
        smlawb(
            fixConst(16.0, 7),
            fixConst(TUNING_NOISE_FLOOR_DB, 7),
            fixConst(0.16, 16)
        )
    )
    var averageGainQ16 = context.averageGainQ16
    val minimumGainQ16 = log2Lin(
        smlawb(
            fixConst(16.0, 7),
            fixConst(TUNING_RELATIVE_MIN_GAIN_DB, 7),
            fixConst(0.16, 16)
        )
    )
    gainOffsetQ16 += smulww(averageGainQ16, minimumGainQ16)

    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        features.gainsQ16[subframeIndex] =
            smulww(features.gainsQ16[subframeIndex], gainMultiplierQ16)
        if (features.gainsQ16[subframeIndex] < 0) {
            features.gainsQ16[subframeIndex] = Int.MAX_VALUE
        }
        features.gainsQ16[subframeIndex] =
            addPosSat32(features.gainsQ16[subframeIndex], gainOffsetQ16)
        averageGainQ16 += smulwb(
            features.gainsQ16[subframeIndex] - averageGainQ16,
            rshiftRound(
                smulbb(speechActivityQ8, fixConst(TUNING_GAIN_SMOOTHING_COEF, 10)),
                2
            )
        )
        subframeIndex++
    }
    context.averageGainQ16 = averageGainQ16
}

private fun applyPreEmphasisGain(
    context: NoiseShapingContext,
    features: FrameFeatures,
    speechActivityQ8: Int
) {
    var preEmphasisGainQ16 = fixConst(1.0, 16) + rshiftRound(
        mla(
            fixConst(TUNING_INPUT_TILT, 26),
            features.codingQualityQ14,
            fixConst(TUNING_HIGH_RATE_INPUT_TILT, 12)
        ),
        10
    )
    if (features.inputTiltQ15 <= 0 && features.signalType == SIG_TYPE_UNVOICED) {
        if (context.samplingRateKHz == 24 || context.samplingRateKHz == 16) {
            val deEssDB = if (context.samplingRateKHz == 24) {
                TUNING_DE_ESSER_COEF_SWB_DB
            } else {
                TUNING_DE_ESSER_COEF_WB_DB
            }
            val deEssStrengthQ15 = smulww(
                -features.inputTiltQ15,
                smulbb(speechActivityQ8, fixConst(1.0, 8) - features.sparsenessQ8)
            )
            val deEssGainQ16 = log2Lin(
                fixConst(16.0, 7) - smulwb(
                    deEssStrengthQ15,
                    smulwb(fixConst(deEssDB, 7), fixConst(0.16, 17))
                )
            )
            preEmphasisGainQ16 = smulww(preEmphasisGainQ16, deEssGainQ16)
        }
    }

    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        features.gainsPreQ14[subframeIndex] =
            smulwb(preEmphasisGainQ16, features.gainsPreQ14[subframeIndex])
        subframeIndex++
    }
}

private fun calculateShapingTilt(
    context: NoiseShapingContext,
    features: FrameFeatures,
    speechActivityQ8: Int
): Int {
    val lowFrequencyStrengthQ16 = multiply(
        fixConst(TUNING_LOW_FREQ_SHAPING, 0),
        fixConst(1.0, 16) + smulbb(
            fixConst(TUNING_LOW_QUALITY_LOW_FREQ_SHAPING_DECR, 1),
            features.inputQualityBandsQ15[0] - fixConst(1.0, 15)
        )
    )
    if (features.signalType == SIG_TYPE_VOICED) {
        val inverseSamplingRateQ14 = div3216(fixConst(0.2, 14), context.samplingRateKHz)
        var subframeIndex = 0
        while (subframeIndex < NB_SUBFR) {
            val lowFrequencyCoefficientQ14 = inverseSamplingRateQ14 +
                    div3216(fixConst(3.0, 14), max(features.pitchL[subframeIndex], 1))
            features.lowFrequencyShapingQ14[subframeIndex] =
                leftShift(
                    fixConst(1.0, 14) - lowFrequencyCoefficientQ14 -
                            smulwb(lowFrequencyStrengthQ16, lowFrequencyCoefficientQ14),
                    16
                ) or
                        ((lowFrequencyCoefficientQ14 - fixConst(1.0, 14)) and 65535)
            subframeIndex++
        }
        return -fixConst(TUNING_HP_NOISE_COEF, 16) - smulwb(
            fixConst(1.0, 16) - fixConst(TUNING_HP_NOISE_COEF, 16),
            smulwb(fixConst(TUNING_HARM_HP_NOISE_COEF, 24), speechActivityQ8)
        )
    }

    val lowFrequencyCoefficientQ14 = div3216(21299, context.samplingRateKHz)
    features.lowFrequencyShapingQ14[0] =
        leftShift(
            fixConst(1.0, 14) - lowFrequencyCoefficientQ14 -
                    smulwb(
                        lowFrequencyStrengthQ16,
                        smulwb(fixConst(0.6, 16), lowFrequencyCoefficientQ14)
                    ),
            16
        ) or
                ((lowFrequencyCoefficientQ14 - fixConst(1.0, 14)) and 65535)
    var subframeIndex = 1
    while (subframeIndex < NB_SUBFR) {
        features.lowFrequencyShapingQ14[subframeIndex] = features.lowFrequencyShapingQ14[0]
        subframeIndex++
    }
    return -fixConst(TUNING_HP_NOISE_COEF, 16)
}

private fun smoothShapingControls(
    shapeState: ShapeState,
    features: FrameFeatures,
    pitchAnalysis: PitchAnalysisContext,
    tiltQ16: Int
) {
    var harmonicBoostQ16 = smulwb(
        smulwb(
            fixConst(1.0, 17) - leftShift(features.codingQualityQ14, 3),
            pitchAnalysis.ltpCorrelationQ15
        ),
        fixConst(TUNING_LOW_RATE_HARMONIC_BOOST, 16)
    )
    harmonicBoostQ16 = smlawb(
        harmonicBoostQ16,
        fixConst(1.0, 16) - leftShift(features.inputQualityQ14, 2),
        fixConst(TUNING_LOW_INPUT_QUALITY_HARMONIC_BOOST, 16)
    )

    val harmonicShapeGainQ16 = if (features.signalType == SIG_TYPE_VOICED) {
        var shapedGainQ16 = smlawb(
            fixConst(TUNING_HARMONIC_SHAPING, 16),
            fixConst(1.0, 16) - smulwb(
                fixConst(1.0, 18) - leftShift(features.codingQualityQ14, 4),
                features.inputQualityQ14
            ),
            fixConst(TUNING_HIGH_RATE_OR_LOW_QUALITY_HARMONIC_SHAPING, 16)
        )
        shapedGainQ16 = smulwb(
            leftShift(shapedGainQ16, 1),
            sqrtApprox(leftShift(pitchAnalysis.ltpCorrelationQ15, 15))
        )
        shapedGainQ16
    } else {
        0
    }

    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        shapeState.harmonicBoostSmoothingQ16 = smlawb(
            shapeState.harmonicBoostSmoothingQ16,
            harmonicBoostQ16 - shapeState.harmonicBoostSmoothingQ16,
            fixConst(TUNING_SUBFR_SMTH_COEF, 16)
        )
        shapeState.harmonicShapeGainSmoothingQ16 = smlawb(
            shapeState.harmonicShapeGainSmoothingQ16,
            harmonicShapeGainQ16 - shapeState.harmonicShapeGainSmoothingQ16,
            fixConst(TUNING_SUBFR_SMTH_COEF, 16)
        )
        shapeState.tiltSmoothingQ16 = smlawb(
            shapeState.tiltSmoothingQ16,
            tiltQ16 - shapeState.tiltSmoothingQ16,
            fixConst(TUNING_SUBFR_SMTH_COEF, 16)
        )
        features.harmBoostQ14[subframeIndex] = rshiftRound(shapeState.harmonicBoostSmoothingQ16, 2)
        features.harmShapeGainQ14[subframeIndex] =
            rshiftRound(shapeState.harmonicShapeGainSmoothingQ16, 2)
        features.tiltQ14[subframeIndex] = rshiftRound(shapeState.tiltSmoothingQ16, 2)
        subframeIndex++
    }
}

private fun analyzeNoiseShapeQuality(
    context: NoiseShapingContext,
    pitchAnalysis: PitchAnalysisContext,
    features: FrameFeatures,
    pitchResidual: IntArray,
    pitchResidualOffset: Int,
    speechActivityQ8: Int,
    targetSnrQ7: Int,
    bufferedChannelMilliseconds: Int,
    inBandFecSnrCompQ8: Int
): Int {
    val buffers = context.buffers
    val energyScale = buffers.energyScale
    val residualEnergy = buffers.residualEnergy

    features.currentSnrDbQ7 = targetSnrQ7 - smulwb(
        leftShift(bufferedChannelMilliseconds, 7),
        fixConst(0.05, 16)
    )
    if (speechActivityQ8 > fixConst(TUNING_LBRR_SPEECH_ACTIVITY_THRESHOLD, 8)) {
        features.currentSnrDbQ7 -= rightShift(inBandFecSnrCompQ8, 1)
    }

    features.inputQualityQ14 = rightShift(
        features.inputQualityBandsQ15[0] + features.inputQualityBandsQ15[1],
        2
    )
    features.codingQualityQ14 = rightShift(
        sigmoidQ15(rshiftRound(features.currentSnrDbQ7 - fixConst(18.0, 7), 4)),
        1
    )

    var backgroundProbabilityQ8 = fixConst(1.0, 8) - speechActivityQ8
    backgroundProbabilityQ8 = smulwb(leftShift(backgroundProbabilityQ8, 8), backgroundProbabilityQ8)
    var snrAdjustmentDBQ7 = smlawb(
        features.currentSnrDbQ7,
        smulbb(fixConst(-TUNING_BG_SNR_DECR_DB, 7) shr 5, backgroundProbabilityQ8),
        smulwb(fixConst(1.0, 14) + features.inputQualityQ14, features.codingQualityQ14)
    )

    if (features.signalType == SIG_TYPE_VOICED) {
        snrAdjustmentDBQ7 = smlawb(
            snrAdjustmentDBQ7,
            fixConst(TUNING_HARM_SNR_INCR_DB, 8),
            pitchAnalysis.ltpCorrelationQ15
        )
        features.quantOffsetType = 0
        features.sparsenessQ8 = 0
    } else {
        snrAdjustmentDBQ7 = smlawb(
            snrAdjustmentDBQ7,
            smlawb(fixConst(6.0, 9), -fixConst(0.4, 18), features.currentSnrDbQ7),
            fixConst(1.0, 14) - features.inputQualityQ14
        )

        val pitchSegmentSampleCount = leftShift(context.samplingRateKHz, 1)
        var energyVariationQ7 = 0
        var previousLogEnergyQ7 = 0
        var pitchResidualIndex = pitchResidualOffset
        var segmentIndex = 0
        while (segmentIndex < FRAME_LENGTH_MS / 2) {
            sumSquaresShift(
                residualEnergy,
                energyScale,
                pitchResidual,
                pitchResidualIndex,
                pitchSegmentSampleCount
            )
            val normalizedEnergy =
                residualEnergy.value + rightShift(pitchSegmentSampleCount, energyScale.value)
            val logEnergyQ7 = lin2Log(max(normalizedEnergy, 1))
            if (segmentIndex > 0) {
                energyVariationQ7 += abs(logEnergyQ7 - previousLogEnergyQ7)
            }
            previousLogEnergyQ7 = logEnergyQ7
            pitchResidualIndex += pitchSegmentSampleCount
            segmentIndex++
        }

        val sparsenessQ8 = rightShift(
            sigmoidQ15(
                smulwb(energyVariationQ7 - fixConst(5.0, 7), fixConst(0.1, 16))
            ),
            7
        )
        features.sparsenessQ8 = sparsenessQ8
        features.quantOffsetType =
            if (sparsenessQ8 > fixConst(TUNING_SPARSENESS_THRESHOLD_QNT_OFFSET, 8)) 0 else 1
        snrAdjustmentDBQ7 = smlawb(
            snrAdjustmentDBQ7,
            fixConst(TUNING_SPARSE_SNR_INCR_DB, 15),
            sparsenessQ8 - fixConst(0.5, 8)
        )
    }
    return snrAdjustmentDBQ7
}

@Suppress("DuplicatedCode")
private fun analyzeShapingCoefficients(
    context: NoiseShapingContext,
    features: FrameFeatures,
    input: IntArray,
    inputOffset: Int
) {
    val buffers = context.buffers
    val energyScale = buffers.energyScale
    val residualEnergy = buffers.residualEnergy
    val autocorrelation = buffers.autocorrelation
    val reflectionCoefficientsQ16 = buffers.reflectionCoefficientsQ16
    val ar1Q24 = buffers.ar1Q24
    val ar2Q24 = buffers.ar2Q24
    val windowedInput = buffers.windowedInput
    var windowOffset = inputOffset - context.lookaheadShape

    val noiseStrengthQ16 = smulwb(
        features.predGainQ16,
        fixConst(TUNING_FIND_PITCH_WHITE_NOISE_FRACTION, 16)
    )
    val baseBandwidthExpansionQ16 = div32VarQ(
        fixConst(TUNING_BANDWIDTH_EXPANSION, 16),
        smlaww(fixConst(1.0, 16), noiseStrengthQ16, noiseStrengthQ16),
        16
    )
    val bandwidthExpansionDeltaQ16 = smulwb(
        fixConst(1.0, 16) - smulbb(3, features.codingQualityQ14),
        fixConst(TUNING_LOW_RATE_BANDWIDTH_EXPANSION_DELTA, 16)
    )
    val secondBandwidthExpansionQ16 =
        baseBandwidthExpansionQ16 + bandwidthExpansionDeltaQ16
    val firstBandwidthExpansionQ16 = div3216(
        leftShift(baseBandwidthExpansionQ16 - bandwidthExpansionDeltaQ16, 14),
        rightShift(secondBandwidthExpansionQ16, 2)
    )
    val warpingQ16 = if (context.warpingQ16 > 0) {
        smlawb(context.warpingQ16, features.codingQualityQ14, fixConst(0.01, 18))
    } else {
        0
    }

    val shapingLpcOrder = context.shapingLpcOrder
    val shapeWindowLength = context.shapeWindowLength
    val flatPartSampleCount = context.samplingRateKHz * 5
    val slopePartSampleCount = rightShift(shapeWindowLength - flatPartSampleCount, 1)
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        applySineWindow(windowedInput, 0, input, windowOffset, 1, slopePartSampleCount)
        var flatSampleIndex = 0
        while (flatSampleIndex < flatPartSampleCount) {
            windowedInput[slopePartSampleCount + flatSampleIndex] =
                input[windowOffset + slopePartSampleCount + flatSampleIndex]
            flatSampleIndex++
        }
        applySineWindow(
            windowedInput,
            slopePartSampleCount + flatPartSampleCount,
            input,
            windowOffset + slopePartSampleCount + flatPartSampleCount,
            2,
            slopePartSampleCount
        )
        windowOffset += context.subframeLength

        if (warpingQ16 > 0) {
            warpedAutocorrelation(
                buffers.warpedFilterStateQ14,
                buffers.warpedCorrelationQ10,
                autocorrelation,
                energyScale,
                windowedInput,
                warpingQ16,
                shapeWindowLength,
                shapingLpcOrder
            )
        } else {
            autocorr(
                autocorrelation,
                energyScale,
                windowedInput,
                0,
                shapeWindowLength,
                shapingLpcOrder + 1
            )
        }

        autocorrelation[0] += max(
            smulwb(
                rightShift(autocorrelation[0], 4),
                fixConst(TUNING_SHAPE_WHITE_NOISE_FRACTION, 20)
            ),
            1
        )

        residualEnergy.value = schur64(reflectionCoefficientsQ16, autocorrelation, shapingLpcOrder)
        ar2Q24.fill(0)
        k2AQ16(ar2Q24, reflectionCoefficientsQ16, shapingLpcOrder)

        var energyScaleExponent = -energyScale.value
        if (energyScaleExponent and 1 != 0) {
            energyScaleExponent -= 1
            residualEnergy.value = rightShift(residualEnergy.value, 1)
        }
        val residualGain = sqrtApprox(max(residualEnergy.value, 1))
        energyScaleExponent = rightShift(energyScaleExponent, 1)
        features.gainsQ16[subframeIndex] = lshiftSat32(residualGain, 16 - energyScaleExponent)

        if (warpingQ16 > 0) {
            val warpedGainQ16 = warpedGain(ar2Q24, warpingQ16, shapingLpcOrder)
            features.gainsQ16[subframeIndex] =
                smulww(features.gainsQ16[subframeIndex], warpedGainQ16)
            if (features.gainsQ16[subframeIndex] < 0) {
                features.gainsQ16[subframeIndex] = Int.MAX_VALUE
            }
        }

        expandBandwidth32(ar2Q24, shapingLpcOrder, secondBandwidthExpansionQ16)
        ar2Q24.copyInto(ar1Q24)
        expandBandwidth32(ar1Q24, shapingLpcOrder, firstBandwidthExpansionQ16)

        val secondInverseGainQ30 = lpcInversePredGainQ24(ar2Q24, shapingLpcOrder)
        val firstInverseGainQ30 = lpcInversePredGainQ24(ar1Q24, shapingLpcOrder)
        val scaledSecondInverseGainQ30 =
            leftShift(smulwb(secondInverseGainQ30, fixConst(0.7, 15)), 1)
        features.gainsPreQ14[subframeIndex] = fixConst(0.3, 14) +
                div32VarQ(scaledSecondInverseGainQ30, max(firstInverseGainQ30, 1), 14)

        limitWarpedCoefficients(
            ar2Q24,
            ar1Q24,
            warpingQ16,
            fixConst(3.999, 24),
            shapingLpcOrder
        )

        var coefficientIndex = 0
        while (coefficientIndex < shapingLpcOrder) {
            val coefficientOffset = subframeIndex * MAX_SHAPE_LPC_ORDER + coefficientIndex
            features.ar1Q13[coefficientOffset] = sat16(rshiftRound(ar1Q24[coefficientIndex], 11))
            features.ar2Q13[coefficientOffset] = sat16(rshiftRound(ar2Q24[coefficientIndex], 11))
            coefficientIndex++
        }
        subframeIndex++
    }
}
