package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.SIG_TYPE_VOICED

// Ported from tsilk.
// Source: hp_variable_cutoff.ts

private const val RADIANS_CONSTANT_Q19 = 1482
private const val LOG2_MINIMUM_FREQUENCY_Q7 = 809
private const val MAXIMUM_FREQUENCY_DELTA_Q7 = 51
private const val TRACKING_SMOOTHING_COEFFICIENT_Q16 = 6554
private const val OUTPUT_SMOOTHING_COEFFICIENT_Q16 = 983
private const val QUALITY_BIAS_Q15 = 19661
private const val RESONANCE_COEFFICIENT_Q9 = 471
private const val ONE_Q28 = 268435456
private const val TWO_Q22 = 8388608
private const val MINIMUM_FREQUENCY_HZ = 80
private const val MAXIMUM_FREQUENCY_HZ = 150

internal fun highPassVariableCutoff(
    context: NoiseShapingContext,
    pitchAnalysis: PitchAnalysisContext,
    features: FrameFeatures,
    output: IntArray,
    input: IntArray,
    speechActivityQ8: Int,
) {
    var pitchFrequencyLogQ7: Int
    var frequencyDeltaQ7: Int
    if (pitchAnalysis.previousSignalType == SIG_TYPE_VOICED) {
        val pitchFrequencyQ16 = div3216(
            leftShift(multiply(context.samplingRateKHz, 1000), 16),
            pitchAnalysis.previousLag,
        )
        pitchFrequencyLogQ7 = lin2Log(pitchFrequencyQ16) - (16 shl 7)
        val qualityQ15 = features.inputQualityBandsQ15[0]
        pitchFrequencyLogQ7 = sub32(
            pitchFrequencyLogQ7,
            smulwb(
                smulwb(leftShift(qualityQ15, 2), qualityQ15),
                pitchFrequencyLogQ7 - LOG2_MINIMUM_FREQUENCY_Q7,
            ),
        )
        pitchFrequencyLogQ7 = add32(
            pitchFrequencyLogQ7,
            rightShift(QUALITY_BIAS_Q15 - qualityQ15, 9),
        )
        frequencyDeltaQ7 = pitchFrequencyLogQ7 - rightShift(context.highPassCutoffTrackingQ15, 8)
        if (frequencyDeltaQ7 < 0) {
            frequencyDeltaQ7 = multiply(frequencyDeltaQ7, 3)
        }
        frequencyDeltaQ7 = limit32(
            frequencyDeltaQ7,
            -MAXIMUM_FREQUENCY_DELTA_Q7,
            MAXIMUM_FREQUENCY_DELTA_Q7,
        )
        context.highPassCutoffTrackingQ15 = smlawb(
            context.highPassCutoffTrackingQ15,
            multiply(leftShift(speechActivityQ8, 1), frequencyDeltaQ7),
            TRACKING_SMOOTHING_COEFFICIENT_Q16,
        )
    }

    context.highPassCutoffSmoothedQ15 = smlawb(
        context.highPassCutoffSmoothedQ15,
        context.highPassCutoffTrackingQ15 - context.highPassCutoffSmoothedQ15,
        OUTPUT_SMOOTHING_COEFFICIENT_Q16,
    )
    features.pitchFreqLowHz = log2Lin(rightShift(context.highPassCutoffSmoothedQ15, 8))
    features.pitchFreqLowHz = limit32(features.pitchFreqLowHz, MINIMUM_FREQUENCY_HZ, MAXIMUM_FREQUENCY_HZ)

    val cutoffQ19 = div3216(
        smulbb(RADIANS_CONSTANT_Q19, features.pitchFreqLowHz),
        context.samplingRateKHz,
    )
    val resonanceQ28 = ONE_Q28 - multiply(RESONANCE_COEFFICIENT_Q9, cutoffQ19)
    val numeratorQ28 = context.highPassNumeratorQ28
    numeratorQ28[0] = resonanceQ28
    numeratorQ28[1] = leftShift(-resonanceQ28, 1)
    numeratorQ28[2] = resonanceQ28
    val resonanceQ22 = rightShift(resonanceQ28, 6)
    val denominatorQ28 = context.highPassDenominatorQ28
    denominatorQ28[0] = smulww(resonanceQ22, smulww(cutoffQ19, cutoffQ19) - TWO_Q22)
    denominatorQ28[1] = smulww(resonanceQ22, resonanceQ22)
    applyBiquadFilter(
        input,
        0,
        numeratorQ28,
        denominatorQ28,
        context.highPassFilterState,
        0,
        output,
        0,
        context.frameLength,
    )
}
