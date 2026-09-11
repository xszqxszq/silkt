package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.util.RefInt

internal const val VAD_N_BANDS: Int = 4

private const val VAD_INTERNAL_SUBFRAMES_LOG2: Int = 2
private const val VAD_INTERNAL_SUBFRAMES: Int = 1 shl VAD_INTERNAL_SUBFRAMES_LOG2
private const val VAD_NOISE_LEVEL_SMOOTH_COEF_Q16: Int = 1024
private const val VAD_NOISE_LEVELS_BIAS: Int = 50
private const val VAD_NEGATIVE_OFFSET_Q5: Int = 128
private const val VAD_SNR_FACTOR_Q16: Int = 45000
private const val VAD_SNR_SMOOTH_COEF_Q18: Int = 4096

internal class VoiceActivityState(
    val analysisFilterState0: IntArray = IntArray(2),
    val analysisFilterState1: IntArray = IntArray(2),
    val analysisFilterState2: IntArray = IntArray(2),
    val lastSubframeEnergies: IntArray = IntArray(VAD_N_BANDS),
    val smoothedEnergyToNoiseRatiosQ8: IntArray = IntArray(VAD_N_BANDS),
    var highPassState: Int = 0,
    val noiseLevels: IntArray = IntArray(VAD_N_BANDS),
    val inverseNoiseLevels: IntArray = IntArray(VAD_N_BANDS),
    val noiseLevelBias: IntArray = IntArray(VAD_N_BANDS),
    var frameCounter: Int = 0,
)

// Ported from tsilk.
// Source: vad.ts

private const val ALL_PASS_EVEN_SECTION_COEFFICIENT = (20623 shl 1)

private const val ALL_PASS_ODD_SECTION_COEFFICIENT = 5394 shl 1

@Suppress("SameParameterValue")
private fun splitIntoBands(
    inData: IntArray,
    inOffset: Int,
    filterState: IntArray,
    outL: IntArray,
    outLOffset: Int,
    outH: IntArray,
    outHOffset: Int,
    sampleCount: Int
) {
    val halfSampleCount = rightShift(sampleCount, 1)
    var sampleIndex = 0
    while (sampleIndex < halfSampleCount) {
        var inputSample = leftShift(inData[inOffset + 2 * sampleIndex], 10)
        val evenSectionDelta = sub32(inputSample, filterState[0])
        val evenSectionCorrection = smlawb(
            evenSectionDelta,
            evenSectionDelta,
            ALL_PASS_EVEN_SECTION_COEFFICIENT
        )
        val evenSectionOutput = add32(
            filterState[0],
            evenSectionCorrection
        )
        filterState[0] = add32(inputSample, evenSectionCorrection)

        inputSample = leftShift(inData[inOffset + 2 * sampleIndex + 1], 10)
        val oddSectionDelta = sub32(inputSample, filterState[1])
        val oddSectionCorrection =
            smulwb(oddSectionDelta, ALL_PASS_ODD_SECTION_COEFFICIENT)
        val oddSectionOutput = add32(
            filterState[1],
            oddSectionCorrection
        )
        filterState[1] = add32(inputSample, oddSectionCorrection)

        outL[outLOffset + sampleIndex] = sat16(
            rshiftRound(add32(oddSectionOutput, evenSectionOutput), 11)
        )
        outH[outHOffset + sampleIndex] = sat16(
            rshiftRound(sub32(oddSectionOutput, evenSectionOutput), 11)
        )
        sampleIndex++
    }
}

internal fun initializeVoiceActivity(state: VoiceActivityState) {
    var band = 0
    while (band < VAD_N_BANDS) {
        state.noiseLevelBias[band] =
            max32(div3216(VAD_NOISE_LEVELS_BIAS, band + 1), 1)
        band++
    }

    band = 0
    while (band < VAD_N_BANDS) {
        state.noiseLevels[band] = multiply(100, state.noiseLevelBias[band])
        state.inverseNoiseLevels[band] = div32(2147483647, state.noiseLevels[band])
        band++
    }

    state.frameCounter = 15
    band = 0
    while (band < VAD_N_BANDS) {
        state.smoothedEnergyToNoiseRatiosQ8[band] = 100 * 256
        band++
    }
}

private val tiltWeights: IntArray = intArrayOf(30000, 6000, -12000, -12000)

private fun updateNoiseLevels(input: IntArray, state: VoiceActivityState) {
    var minimumCoefficient = 0
    if (state.frameCounter < 1000) {
        minimumCoefficient = div3216(32767, rightShift(state.frameCounter, 4) + 1)
    }

    var band = 0
    while (band < VAD_N_BANDS) {
        var noiseLevel = state.noiseLevels[band]
        val energy = addPosSat32(input[band], state.noiseLevelBias[band])
        val inverseEnergy = div32(2147483647, energy)

        val coefficient = when {
            energy > leftShift(noiseLevel, 3) ->
                VAD_NOISE_LEVEL_SMOOTH_COEF_Q16 shr 3

            energy < noiseLevel -> VAD_NOISE_LEVEL_SMOOTH_COEF_Q16
            else -> smulwb(
                smulww(inverseEnergy, noiseLevel),
                VAD_NOISE_LEVEL_SMOOTH_COEF_Q16 shl 1
            )
        }

        state.inverseNoiseLevels[band] = smlawb(
            state.inverseNoiseLevels[band],
            inverseEnergy - state.inverseNoiseLevels[band],
            maxInt(coefficient, minimumCoefficient)
        )
        noiseLevel = minInt(div32(2147483647, state.inverseNoiseLevels[band]), 16777215)
        state.noiseLevels[band] = noiseLevel
        band++
    }
    state.frameCounter++
}

@Suppress("SameReturnValue")
internal fun getSpeechActivityQ8(
    state: VoiceActivityState,
    workspace: VoiceActivityWorkspace,
    speechActivityQ8: RefInt,
    snrDbQ7: RefInt,
    qualityQ15: IntArray,
    tiltQ15: RefInt,
    input: IntArray,
    frameLength: Int
): Int {
    val bands = workspace.bands
    val bandEnergies = workspace.bandEnergies
    val energyToNoiseRatiosQ8 = workspace.energyToNoiseRatiosQ8

    splitIntoBands(input, 0, state.analysisFilterState0, bands[0], 0, bands[3], 0, frameLength)
    splitIntoBands(
        bands[0], 0, state.analysisFilterState1, bands[0], 0, bands[2], 0,
        rightShift(frameLength, 1)
    )
    splitIntoBands(
        bands[0], 0, state.analysisFilterState2, bands[0], 0, bands[1], 0,
        rightShift(frameLength, 2)
    )
    var decimatedFrameLength = rightShift(frameLength, 3)
    bands[0][decimatedFrameLength - 1] = rightShift(bands[0][decimatedFrameLength - 1], 1)
    val highPassState = bands[0][decimatedFrameLength - 1]
    var sampleIndex = decimatedFrameLength - 1
    while (sampleIndex > 0) {
        bands[0][sampleIndex - 1] = rightShift(bands[0][sampleIndex - 1], 1)
        bands[0][sampleIndex] -= bands[0][sampleIndex - 1]
        sampleIndex--
    }
    bands[0][0] -= state.highPassState
    state.highPassState = highPassState

    var band = 0
    while (band < VAD_N_BANDS) {
        decimatedFrameLength = rightShift(
            frameLength,
            minInt(VAD_N_BANDS - band, VAD_N_BANDS - 1)
        )
        val decimatedSubframeLength = rightShift(decimatedFrameLength, VAD_INTERNAL_SUBFRAMES_LOG2)
        var decimatedOffset = 0
        bandEnergies[band] = state.lastSubframeEnergies[band]
        var subframeEnergy = 0

        var subframeIndex = 0
        while (subframeIndex < VAD_INTERNAL_SUBFRAMES) {
            subframeEnergy = 0
            var sampleIndexInSubframe = 0
            while (sampleIndexInSubframe < decimatedSubframeLength) {
                val sample = rightShift(
                    bands[band][sampleIndexInSubframe + decimatedOffset],
                    3
                )
                subframeEnergy = smlabb(subframeEnergy, sample, sample)
                sampleIndexInSubframe++
            }
            if (subframeIndex < VAD_INTERNAL_SUBFRAMES - 1) {
                bandEnergies[band] = addPosSat32(bandEnergies[band], subframeEnergy)
            } else {
                bandEnergies[band] =
                    addPosSat32(bandEnergies[band], rightShift(subframeEnergy, 1))
            }
            decimatedOffset += decimatedSubframeLength
            subframeIndex++
        }
        state.lastSubframeEnergies[band] = subframeEnergy
        band++
    }

    updateNoiseLevels(bandEnergies, state)

    var squaredSnrSum = 0
    var inputTilt = 0
    band = 0
    while (band < VAD_N_BANDS) {
        val speechEnergy = bandEnergies[band] - state.noiseLevels[band]
        if (speechEnergy > 0) {
            if ((bandEnergies[band] and 4286578688.toInt()) == 0) {
                energyToNoiseRatiosQ8[band] = div32(
                    leftShift(bandEnergies[band], 8),
                    state.noiseLevels[band] + 1
                )
            } else {
                energyToNoiseRatiosQ8[band] = div32(
                    bandEnergies[band],
                    rightShift(state.noiseLevels[band], 8) + 1
                )
            }

            var snrQ7 = lin2Log(energyToNoiseRatiosQ8[band]) - 8 * 128
            squaredSnrSum = smlabb(squaredSnrSum, snrQ7, snrQ7)
            if (speechEnergy < 1 shl 20) {
                snrQ7 = smulwb(leftShift(sqrtApprox(speechEnergy), 6), snrQ7)
            }
            inputTilt = smlawb(inputTilt, tiltWeights[band], snrQ7)
        } else {
            energyToNoiseRatiosQ8[band] = 256
        }
        band++
    }

    squaredSnrSum = div3216(squaredSnrSum, VAD_N_BANDS)
    snrDbQ7.value = 3 * sqrtApprox(squaredSnrSum)

    var speechActivityQ15 = sigmoidQ15(
        smulwb(VAD_SNR_FACTOR_Q16, snrDbQ7.value) - VAD_NEGATIVE_OFFSET_Q5
    )
    tiltQ15.value = leftShift(sigmoidQ15(inputTilt) - 16384, 1)

    var weightedSpeechEnergy = 0
    band = 0
    while (band < VAD_N_BANDS) {
        weightedSpeechEnergy +=
            (band + 1) * rightShift(bandEnergies[band] - state.noiseLevels[band], 4)
        band++
    }
    if (weightedSpeechEnergy <= 0) {
        speechActivityQ15 = rightShift(speechActivityQ15, 1)
    } else if (weightedSpeechEnergy < 32768) {
        weightedSpeechEnergy = sqrtApprox(leftShift(weightedSpeechEnergy, 15))
        speechActivityQ15 = smulwb(32768 + weightedSpeechEnergy, speechActivityQ15)
    }

    speechActivityQ8.value = minInt(rightShift(speechActivityQ15, 7), 255)
    val smoothingCoefficientQ16 =
        smulwb(VAD_SNR_SMOOTH_COEF_Q18, smulwb(speechActivityQ15, speechActivityQ15))

    band = 0
    while (band < VAD_N_BANDS) {
        state.smoothedEnergyToNoiseRatiosQ8[band] = smlawb(
            state.smoothedEnergyToNoiseRatiosQ8[band],
            energyToNoiseRatiosQ8[band] - state.smoothedEnergyToNoiseRatiosQ8[band],
            smoothingCoefficientQ16
        )
        val snrQ7 = 3 * (lin2Log(state.smoothedEnergyToNoiseRatiosQ8[band]) - 8 * 128)
        qualityQ15[band] = sigmoidQ15(rightShift(snrQ7 - 16 * 128, 4))
        band++
    }

    return 0
}
