package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.*
import xyz.xszq.silkt.internal.signal.FrameFeatures
import xyz.xszq.silkt.internal.signal.movingAveragePrediction
import xyz.xszq.silkt.internal.tables.ExcitationTables.quantizationOffsetsQ10
import xyz.xszq.silkt.internal.util.RefInt

private fun clampPulse(pulse: Int): Int {
    return pulse.coerceIn(-64, 64)
}

@Suppress("SameParameterValue")
private fun addLShift32(a: Int, b: Int, shift: Int): Int = add32(a, leftShift(b, shift))

@Suppress("SameParameterValue")
private fun subLShift32(a: Int, b: Int, shift: Int): Int = sub32(a, leftShift(b, shift))

@Suppress("SameParameterValue")
private fun subRShift32(a: Int, b: Int, shift: Int): Int = sub32(a, rightShift(b, shift))

private fun copySampleState(destination: SampleState, source: SampleState) {
    destination.pulseQ10 = source.pulseQ10
    destination.rateDistortionQ10 = source.rateDistortionQ10
    destination.quantizedSignalQ14 = source.quantizedSignalQ14
    destination.lowFrequencyShapingStateQ12 = source.lowFrequencyShapingStateQ12
    destination.ltpShapeQ10 = source.ltpShapeQ10
    destination.lpcExcitationQ16 = source.lpcExcitationQ16
}

private fun selectMinimumDistortionIndex(
    delayedDecisions: Array<DelayedDecisionState>,
    delayedDecisionCount: Int
): Int {
    var winnerIndex = 0
    var minimumRateDistortionQ10 = delayedDecisions[0].rateDistortionQ10
    var searchIndex = 1
    while (searchIndex < delayedDecisionCount) {
        val rateDistortionQ10 = delayedDecisions[searchIndex].rateDistortionQ10
        if (rateDistortionQ10 < minimumRateDistortionQ10) {
            minimumRateDistortionQ10 = rateDistortionQ10
            winnerIndex = searchIndex
        }
        searchIndex++
    }
    return winnerIndex
}

private fun penalizeRejectedDecisions(
    delayedDecisions: Array<DelayedDecisionState>,
    delayedDecisionCount: Int,
    winnerIndex: Int
) {
    var decisionIndex = 0
    while (decisionIndex < delayedDecisionCount) {
        if (decisionIndex != winnerIndex) {
            delayedDecisions[decisionIndex].rateDistortionQ10 =
                add32(delayedDecisions[decisionIndex].rateDistortionQ10, INT32_MAX shr 4)
        }
        decisionIndex++
    }
}

private fun flushDelayedDecisionSamples(
    state: NoiseShapeQuantizerState,
    winner: DelayedDecisionState,
    quantizedPulses: IntArray,
    pulseOffset: Int,
    quantizedSignalOffset: Int,
    sampleHistoryIndex: Int,
    decisionDelay: Int
) {
    var lastSampleIndex = sampleHistoryIndex + decisionDelay
    var delayedSampleIndex = 0
    while (delayedSampleIndex < decisionDelay) {
        lastSampleIndex = (lastSampleIndex - 1) and DECISION_DELAY_MASK
        quantizedPulses[pulseOffset + delayedSampleIndex - decisionDelay] =
            rightShift(winner.pulseQ10[lastSampleIndex], 10)
        state.quantizedSignalQ10[quantizedSignalOffset + delayedSampleIndex - decisionDelay] =
            sat16(
                rshiftRound(
                    smulww(winner.quantizedSignalQ10[lastSampleIndex], winner.gainQ16[lastSampleIndex]),
                    10
                )
            )
        state.ltpShapeQ10[state.ltpShapeWriteIndex - decisionDelay + delayedSampleIndex] =
            winner.ltpShapeQ10[lastSampleIndex]
        delayedSampleIndex++
    }
}

private fun rewhiteningSubframeMask(interpolationFlag: Int): Int = 3 - (interpolationFlag shl 1)

private fun packHarmonicShapeFirQ14(harmonicShapeGainQ14: Int): Int {
    val lowPartQ14 = rightShift(harmonicShapeGainQ14, 2)
    val highPartQ14 = leftShift(rightShift(harmonicShapeGainQ14, 1), 16)
    return lowPartQ14 or highPartQ14
}

private fun prepareRewhitening(
    state: NoiseShapeQuantizerState,
    frameLength: Int,
    lag: Int,
    predictionLpcOrder: Int,
    predCoefQ12: IntArray,
    predictionCoefficientOffset: Int,
    lpcFilterState: IntArray,
    ltpExcitation: IntArray,
    quantizedSignalStartIndex: Int
) {
    val historyStartIndex = frameLength - lag - predictionLpcOrder - (LTP_ORDER shr 1)
    lpcFilterState.fill(0, 0, predictionLpcOrder)
    movingAveragePrediction(
        state.quantizedSignalQ10,
        quantizedSignalStartIndex,
        predCoefQ12,
        predictionCoefficientOffset,
        lpcFilterState,
        0,
        ltpExcitation,
        historyStartIndex,
        frameLength - historyStartIndex,
        predictionLpcOrder
    )
    state.ltpExcitationWriteIndex = frameLength
    state.shouldRewhiten = true
}

private fun calculateEffectiveDecisionDelay(
    features: FrameFeatures,
    previousLag: Int,
    subframeLength: Int
): Int {
    var decisionDelay = min(DECISION_DELAY, subframeLength)
    if (features.signalType == SIG_TYPE_VOICED) {
        var pitchSubframeIndex = 0
        while (pitchSubframeIndex < NB_SUBFR) {
            decisionDelay = min(decisionDelay, features.pitchL[pitchSubframeIndex] - (LTP_ORDER shr 1) - 1)
            pitchSubframeIndex++
        }
    } else if (previousLag > 0) {
        decisionDelay = min(decisionDelay, previousLag - (LTP_ORDER shr 1) - 1)
    }
    return decisionDelay
}

private fun flushDelayedDecisionExcitation(
    state: NoiseShapeQuantizerState,
    winner: DelayedDecisionState,
    scaledLtpExcitationQ16: IntArray,
    sampleHistoryIndex: Int,
    decisionDelay: Int
) {
    var excitationFlushIndex = 0
    while (excitationFlushIndex < decisionDelay) {
        val winnerSampleIndex =
            (sampleHistoryIndex + decisionDelay - excitationFlushIndex - 1) and DECISION_DELAY_MASK
        val outputIndex = state.ltpExcitationWriteIndex - decisionDelay + excitationFlushIndex
        scaledLtpExcitationQ16[outputIndex] = winner.lpcExcitationQ16[winnerSampleIndex]
        excitationFlushIndex++
    }
}

private fun copyWinnerQuantizerState(
    state: NoiseShapeQuantizerState,
    winner: DelayedDecisionState,
    subframeLength: Int
) {
    winner.lpcStateQ14.copyInto(
        state.lpcStateQ14,
        destinationOffset = 0,
        startIndex = subframeLength,
        endIndex = subframeLength + LPC_HISTORY_BUFFER_LENGTH
    )
    winner.autoregressiveStateQ14.copyInto(
        state.autoregressiveStateQ14,
        destinationOffset = 0,
        startIndex = 0,
        endIndex = MAX_SHAPE_LPC_ORDER
    )
    state.lowFrequencyShapingStateQ12 = winner.lowFrequencyShapingStateQ12
}

private fun initializeDelayedDecisionStates(
    state: NoiseShapeQuantizerState,
    delayedDecisions: Array<DelayedDecisionState>,
    delayedDecisionCount: Int,
    frameLength: Int,
    frameStartSeed: Int
) {
    var decisionIndex = 0
    while (decisionIndex < delayedDecisionCount) {
        val decision = delayedDecisions[decisionIndex]
        decision.seed = (decisionIndex + frameStartSeed) and 3
        decision.frameStartSeed = decision.seed
        decision.rateDistortionQ10 = 0
        decision.lowFrequencyShapingStateQ12 = state.lowFrequencyShapingStateQ12
        decision.seedHistory.fill(0)
        decision.pulseQ10.fill(0)
        decision.quantizedSignalQ10.fill(0)
        decision.lpcExcitationQ16.fill(0)
        decision.ltpShapeQ10.fill(0)
        decision.gainQ16.fill(0)
        decision.ltpShapeQ10[0] = state.ltpShapeQ10[frameLength - 1]

        var sourceIndex = 0
        while (sourceIndex < LPC_HISTORY_BUFFER_LENGTH) {
            decision.lpcStateQ14[sourceIndex] = state.lpcStateQ14[sourceIndex]
            sourceIndex++
        }
        var autoregressiveCopyIndex = 0
        while (autoregressiveCopyIndex < MAX_SHAPE_LPC_ORDER) {
            decision.autoregressiveStateQ14[autoregressiveCopyIndex] =
                state.autoregressiveStateQ14[autoregressiveCopyIndex]
            autoregressiveCopyIndex++
        }
        decisionIndex++
    }
}

private fun selectBestSampleStateIndex(
    sampleStates: Array<Array<SampleState>>,
    delayedDecisionCount: Int
): Int {
    var minimumRateDistortionQ10 = sampleStates[0][0].rateDistortionQ10
    var winnerIndex = 0
    var winnerSearchIndex = 1
    while (winnerSearchIndex < delayedDecisionCount) {
        if (sampleStates[winnerSearchIndex][0].rateDistortionQ10 < minimumRateDistortionQ10) {
            minimumRateDistortionQ10 = sampleStates[winnerSearchIndex][0].rateDistortionQ10
            winnerIndex = winnerSearchIndex
        }
        winnerSearchIndex++
    }
    return winnerIndex
}

private fun penalizeMismatchedSeedStates(
    delayedDecisions: Array<DelayedDecisionState>,
    sampleStates: Array<Array<SampleState>>,
    delayedDecisionCount: Int,
    winnerIndex: Int,
    lastSampleIndex: Int
) {
    val winnerRandState = delayedDecisions[winnerIndex].seedHistory[lastSampleIndex]
    var seedSearchIndex = 0
    while (seedSearchIndex < delayedDecisionCount) {
        if (delayedDecisions[seedSearchIndex].seedHistory[lastSampleIndex] != winnerRandState) {
            sampleStates[seedSearchIndex][0].rateDistortionQ10 =
                add32(sampleStates[seedSearchIndex][0].rateDistortionQ10, INT32_MAX shr 4)
            sampleStates[seedSearchIndex][1].rateDistortionQ10 =
                add32(sampleStates[seedSearchIndex][1].rateDistortionQ10, INT32_MAX shr 4)
        }
        seedSearchIndex++
    }
}

private fun replaceWeakDelayedDecision(
    delayedDecisions: Array<DelayedDecisionState>,
    sampleStates: Array<Array<SampleState>>,
    delayedDecisionCount: Int,
    sampleIndex: Int
) {
    var maximumRateDistortionQ10 = sampleStates[0][0].rateDistortionQ10
    var minimumRateDistortionQ10 = sampleStates[0][1].rateDistortionQ10
    var maximumDistortionIndex = 0
    var minimumDistortionIndex = 0
    var replacementSearchIndex = 1
    while (replacementSearchIndex < delayedDecisionCount) {
        if (sampleStates[replacementSearchIndex][0].rateDistortionQ10 > maximumRateDistortionQ10) {
            maximumRateDistortionQ10 = sampleStates[replacementSearchIndex][0].rateDistortionQ10
            maximumDistortionIndex = replacementSearchIndex
        }
        if (sampleStates[replacementSearchIndex][1].rateDistortionQ10 < minimumRateDistortionQ10) {
            minimumRateDistortionQ10 = sampleStates[replacementSearchIndex][1].rateDistortionQ10
            minimumDistortionIndex = replacementSearchIndex
        }
        replacementSearchIndex++
    }
    if (minimumRateDistortionQ10 < maximumRateDistortionQ10) {
        copyDelayedDecisionState(
            delayedDecisions[maximumDistortionIndex],
            delayedDecisions[minimumDistortionIndex],
            sampleIndex
        )
        copySampleState(
            sampleStates[maximumDistortionIndex][0],
            sampleStates[minimumDistortionIndex][1]
        )
    }
}

private fun commitDelayedDecisionSample(
    decision: DelayedDecisionState,
    selectedState: SampleState,
    sampleIndex: Int,
    historyIndex: Int,
    gainQ16: Int
) {
    decision.lowFrequencyShapingStateQ12 = selectedState.lowFrequencyShapingStateQ12
    decision.lpcStateQ14[LPC_HISTORY_BUFFER_LENGTH + sampleIndex] = selectedState.quantizedSignalQ14
    decision.quantizedSignalQ10[historyIndex] = rightShift(selectedState.quantizedSignalQ14, 4)
    decision.pulseQ10[historyIndex] = selectedState.pulseQ10
    decision.lpcExcitationQ16[historyIndex] = selectedState.lpcExcitationQ16
    decision.ltpShapeQ10[historyIndex] = selectedState.ltpShapeQ10
    decision.seed = addRshift(decision.seed, selectedState.pulseQ10, 10)
    decision.seedHistory[historyIndex] = decision.seed
    decision.rateDistortionQ10 = selectedState.rateDistortionQ10
    decision.gainQ16[historyIndex] = gainQ16
}

private fun resetDelayedDecisionHistories(
    delayedDecisions: Array<DelayedDecisionState>,
    delayedDecisionCount: Int,
    subframeLength: Int
) {
    var decisionResetIndex = 0
    while (decisionResetIndex < delayedDecisionCount) {
        val decision = delayedDecisions[decisionResetIndex]
        var historyResetIndex = 0
        while (historyResetIndex < LPC_HISTORY_BUFFER_LENGTH) {
            decision.lpcStateQ14[historyResetIndex] = decision.lpcStateQ14[subframeLength + historyResetIndex]
            historyResetIndex++
        }
        decisionResetIndex++
    }
}

@Suppress("DuplicatedCode")
private fun quantizeNoiseShapeSubframe(
    state: NoiseShapeQuantizerState,
    signalType: Int,
    scaledInputQ10: IntArray,
    quantizedPulses: IntArray,
    pulseOffset: Int,
    quantizedSignalOutputQ10: IntArray,
    quantizedSignalOffset: Int,
    scaledLtpExcitationQ16: IntArray,
    predictionCoefficientsQ12: IntArray,
    predictionCoefficientOffset: Int,
    ltpCoefficientsQ14: IntArray,
    ltpCoefficientOffset: Int,
    autoregressiveShapingCoefficientsQ13: IntArray,
    autoregressiveShapingCoefficientOffset: Int,
    lag: Int,
    harmShapeFirPackedQ14: Int,
    tiltQ14: Int,
    lowFrequencyShapingQ14: Int,
    gainQ16: Int,
    lambdaQ10: Int,
    offsetQ10: Int,
    subframeLength: Int,
    shapingLpcOrder: Int,
    predictionLpcOrder: Int
) {
    var ltpShapeLagIndex = state.ltpShapeWriteIndex - lag + (HARMONIC_SHAPE_FIR_TAPS shr 1)
    var ltpExcitationLagIndex = state.ltpExcitationWriteIndex - lag + (LTP_ORDER shr 1)
    var lpcHistoryIndex = LPC_HISTORY_BUFFER_LENGTH - 1
    val firstThresholdQ10 = sub32(-1536, rightShift(lambdaQ10, 1))
    var secondThresholdQ10 = sub32(-512, rightShift(lambdaQ10, 1))
    secondThresholdQ10 = addRshift(secondThresholdQ10, smulbb(offsetQ10, lambdaQ10), 10)
    val thirdThresholdQ10 = add32(512, rightShift(lambdaQ10, 1))
    var sampleIndex = 0
    while (sampleIndex < subframeLength) {
        state.randomSeed = rand(state.randomSeed)
        val dither = rightShift(state.randomSeed, 31)
        var lpcPredQ10 = smulwb(
            state.lpcStateQ14[lpcHistoryIndex],
            predictionCoefficientsQ12[predictionCoefficientOffset]
        )
        var predictionCoefficientIndex = 1
        while (predictionCoefficientIndex < predictionLpcOrder) {
            lpcPredQ10 = smlawb(
                lpcPredQ10,
                state.lpcStateQ14[(lpcHistoryIndex - predictionCoefficientIndex)],
                predictionCoefficientsQ12[(predictionCoefficientOffset + predictionCoefficientIndex)]
            )
            predictionCoefficientIndex++
        }
        var ltpPredQ14 = 0
        if (signalType == SIG_TYPE_VOICED) {
            ltpPredQ14 = smulwb(scaledLtpExcitationQ16[ltpExcitationLagIndex], ltpCoefficientsQ14[ltpCoefficientOffset])
            var ltpCoefficientIndex = 1
            while (ltpCoefficientIndex < LTP_ORDER) {
                ltpPredQ14 = smlawb(
                    ltpPredQ14,
                    scaledLtpExcitationQ16[(ltpExcitationLagIndex - ltpCoefficientIndex)],
                    ltpCoefficientsQ14[(ltpCoefficientOffset + ltpCoefficientIndex)]
                )
                ltpCoefficientIndex++
            }
            ltpExcitationLagIndex++
        }
        var currentStateQ14 = state.lpcStateQ14[lpcHistoryIndex]
        var previousStateQ14 = state.autoregressiveStateQ14[0]
        state.autoregressiveStateQ14[0] = currentStateQ14
        var autoregressiveNoiseQ10 = smulwb(
            currentStateQ14,
            autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset]
        )
        var autoregressiveCoefficientIndex = 2
        while (autoregressiveCoefficientIndex < shapingLpcOrder) {
            currentStateQ14 = state.autoregressiveStateQ14[autoregressiveCoefficientIndex - 1]
            state.autoregressiveStateQ14[autoregressiveCoefficientIndex - 1] = previousStateQ14
            autoregressiveNoiseQ10 = smlawb(
                autoregressiveNoiseQ10,
                previousStateQ14,
                autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset + autoregressiveCoefficientIndex - 1]
            )
            previousStateQ14 = state.autoregressiveStateQ14[autoregressiveCoefficientIndex]
            state.autoregressiveStateQ14[autoregressiveCoefficientIndex] = currentStateQ14
            autoregressiveNoiseQ10 = smlawb(
                autoregressiveNoiseQ10,
                currentStateQ14,
                autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset + autoregressiveCoefficientIndex]
            )
            autoregressiveCoefficientIndex += 2
        }
        state.autoregressiveStateQ14[shapingLpcOrder - 1] = previousStateQ14
        autoregressiveNoiseQ10 = smlawb(
            autoregressiveNoiseQ10,
            previousStateQ14,
            autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset + shapingLpcOrder - 1]
        )
        autoregressiveNoiseQ10 = rightShift(autoregressiveNoiseQ10, 1)
        autoregressiveNoiseQ10 = smlawb(autoregressiveNoiseQ10, state.lowFrequencyShapingStateQ12, tiltQ14)
        var lowFrequencyNoiseQ10 = leftShift(
            smulwb(state.ltpShapeQ10[state.ltpShapeWriteIndex - 1], lowFrequencyShapingQ14),
            2
        )
        lowFrequencyNoiseQ10 = smlawt(lowFrequencyNoiseQ10, state.lowFrequencyShapingStateQ12, lowFrequencyShapingQ14)
        var longTermPredictionNoiseQ14 = 0
        if (lag > 0) {
            longTermPredictionNoiseQ14 = smulwb(
                add32(state.ltpShapeQ10[ltpShapeLagIndex], state.ltpShapeQ10[(ltpShapeLagIndex - 2)]),
                harmShapeFirPackedQ14
            )
            longTermPredictionNoiseQ14 =
                smlawt(longTermPredictionNoiseQ14, state.ltpShapeQ10[(ltpShapeLagIndex - 1)], harmShapeFirPackedQ14)
            longTermPredictionNoiseQ14 = leftShift(longTermPredictionNoiseQ14, 6)
            ltpShapeLagIndex++
        }
        var predictionQ10 = sub32(ltpPredQ14, longTermPredictionNoiseQ14)
        predictionQ10 = rightShift(predictionQ10, 4)
        predictionQ10 = add32(predictionQ10, lpcPredQ10)
        predictionQ10 = sub32(predictionQ10, autoregressiveNoiseQ10)
        predictionQ10 = sub32(predictionQ10, lowFrequencyNoiseQ10)
        var residualQ10 = sub32(scaledInputQ10[sampleIndex], predictionQ10)
        residualQ10 = sub32((((residualQ10 xor dither)) - dither), offsetQ10)
        residualQ10 = limit32(residualQ10, (-((64 shl 10))), (64 shl 10))
        var pulseQ0 = 0
        var pulseQ10 = 0
        if (residualQ10 < secondThresholdQ10) {
            if (residualQ10 < firstThresholdQ10) {
                pulseQ0 = rshiftRound(add32(residualQ10, rightShift(lambdaQ10, 1)), 10)
                pulseQ10 = leftShift(pulseQ0, 10)
            } else {
                pulseQ0 = -1
                pulseQ10 = -1024
            }
        } else if (residualQ10 > thirdThresholdQ10) {
            pulseQ0 = rshiftRound(sub32(residualQ10, rightShift(lambdaQ10, 1)), 10)
            pulseQ10 = leftShift(pulseQ0, 10)
        }
        pulseQ0 = clampPulse(pulseQ0)
        quantizedPulses[(pulseOffset + sampleIndex)] = pulseQ0
        var excitationQ10 = add32(pulseQ10, offsetQ10)
        excitationQ10 = (excitationQ10 xor dither) - dither
        val lpcExcitationQ10 = add32(excitationQ10, rshiftRound(ltpPredQ14, 4))
        val quantizedSignalQ10 = add32(lpcExcitationQ10, lpcPredQ10)
        quantizedSignalOutputQ10[(quantizedSignalOffset + sampleIndex)] =
            sat16(rshiftRound(smulww(quantizedSignalQ10, gainQ16), 10))
        lpcHistoryIndex++
        state.lpcStateQ14[lpcHistoryIndex] = leftShift(quantizedSignalQ10, 4)
        val lowFrequencyShapingQ10 = sub32(quantizedSignalQ10, autoregressiveNoiseQ10)
        state.lowFrequencyShapingStateQ12 = leftShift(lowFrequencyShapingQ10, 2)
        state.ltpShapeQ10[state.ltpShapeWriteIndex] = sub32(lowFrequencyShapingQ10, lowFrequencyNoiseQ10)
        scaledLtpExcitationQ16[state.ltpExcitationWriteIndex] = leftShift(lpcExcitationQ10, 6)
        state.ltpShapeWriteIndex++
        state.ltpExcitationWriteIndex++
        state.randomSeed = add32(state.randomSeed, quantizedPulses[(pulseOffset + sampleIndex)])
        sampleIndex++
    }
    var historyIndex = 0
    while (historyIndex < LPC_HISTORY_BUFFER_LENGTH) {
        state.lpcStateQ14[historyIndex] = state.lpcStateQ14[subframeLength + historyIndex]
        historyIndex++
    }
}

@Suppress("DuplicatedCode")
private fun quantizeNoiseShapeSubframeDelayedDecision(
    state: NoiseShapeQuantizerState,
    delayedDecisions: Array<DelayedDecisionState>,
    sampleStates: Array<Array<SampleState>>,
    signalType: Int,
    scaledInputQ10: IntArray,
    quantizedPulses: IntArray,
    pulseOffset: Int,
    quantizedSignalOutputQ10: IntArray,
    quantizedSignalOffset: Int,
    scaledLtpExcitationQ16: IntArray,
    predictionCoefficientsQ12: IntArray,
    predictionCoefficientOffset: Int,
    ltpCoefficientsQ14: IntArray,
    ltpCoefficientOffset: Int,
    autoregressiveShapingCoefficientsQ13: IntArray,
    autoregressiveShapingCoefficientOffset: Int,
    lag: Int,
    harmShapeFirPackedQ14: Int,
    tiltQ14: Int,
    lowFrequencyShapingQ14: Int,
    gainQ16: Int,
    lambdaQ10: Int,
    offsetQ10: Int,
    subframeLength: Int,
    subframeIndex: Int,
    shapingLpcOrder: Int,
    predictionLpcOrder: Int,
    warpingQ16: Int,
    delayedDecisionCount: Int,
    sampleHistoryIndexRef: RefInt,
    decisionDelay: Int
) {
    val stateLtpShapeQ10 = state.ltpShapeQ10
    var ltpShapeWriteIndex = state.ltpShapeWriteIndex
    var ltpExcitationWriteIndex = state.ltpExcitationWriteIndex
    var ltpShapeLagIndex = ltpShapeWriteIndex - lag + (HARMONIC_SHAPE_FIR_TAPS shr 1)
    var ltpExcitationLagIndex = ltpExcitationWriteIndex - lag + (LTP_ORDER shr 1)
    var sampleHistoryIndex = sampleHistoryIndexRef.value
    var sampleIndex = 0
    while (sampleIndex < subframeLength) {
        var ltpPredQ14 = 0
        if (signalType == SIG_TYPE_VOICED) {
            ltpPredQ14 = smulwb(scaledLtpExcitationQ16[ltpExcitationLagIndex], ltpCoefficientsQ14[ltpCoefficientOffset])
            ltpPredQ14 = smlawb(
                ltpPredQ14,
                scaledLtpExcitationQ16[(ltpExcitationLagIndex - 1)],
                ltpCoefficientsQ14[(ltpCoefficientOffset + 1)]
            )
            ltpPredQ14 = smlawb(
                ltpPredQ14,
                scaledLtpExcitationQ16[(ltpExcitationLagIndex - 2)],
                ltpCoefficientsQ14[(ltpCoefficientOffset + 2)]
            )
            ltpPredQ14 = smlawb(
                ltpPredQ14,
                scaledLtpExcitationQ16[(ltpExcitationLagIndex - 3)],
                ltpCoefficientsQ14[(ltpCoefficientOffset + 3)]
            )
            ltpPredQ14 = smlawb(
                ltpPredQ14,
                scaledLtpExcitationQ16[(ltpExcitationLagIndex - 4)],
                ltpCoefficientsQ14[(ltpCoefficientOffset + 4)]
            )
            ltpExcitationLagIndex++
        }
        var longTermPredictionNoiseQ14 = 0
        if (lag > 0) {
            longTermPredictionNoiseQ14 = smulwb(
                add32(
                    stateLtpShapeQ10[ltpShapeLagIndex],
                    stateLtpShapeQ10[(ltpShapeLagIndex - 2)]
                ),
                harmShapeFirPackedQ14
            )
            longTermPredictionNoiseQ14 =
                smlawt(
                    longTermPredictionNoiseQ14,
                    stateLtpShapeQ10[(ltpShapeLagIndex - 1)],
                    harmShapeFirPackedQ14
                )
            longTermPredictionNoiseQ14 = leftShift(longTermPredictionNoiseQ14, 6)
            ltpShapeLagIndex++
        }
        var decisionIndex = 0
        val lpcHistoryIndex = LPC_HISTORY_BUFFER_LENGTH - 1 + sampleIndex
        while (decisionIndex < delayedDecisionCount) {
            val decision = delayedDecisions[decisionIndex]
            val bestState = sampleStates[decisionIndex][0]
            val secondState = sampleStates[decisionIndex][1]
            decision.seed = rand(decision.seed)
            val dither = rightShift(decision.seed, 31)
            val decisionLpcStateQ14 = decision.lpcStateQ14
            var lpcPredQ10 = smulwb(
                decisionLpcStateQ14[lpcHistoryIndex],
                predictionCoefficientsQ12[predictionCoefficientOffset]
            )
            var predictionCoefficientIndex = 1
            while (predictionCoefficientIndex < predictionLpcOrder) {
                lpcPredQ10 = smlawb(
                    lpcPredQ10,
                    decisionLpcStateQ14[(lpcHistoryIndex - predictionCoefficientIndex)],
                    predictionCoefficientsQ12[(predictionCoefficientOffset + predictionCoefficientIndex)]
                )
                predictionCoefficientIndex++
            }
            var currentStateQ14 = smlawb(
                decisionLpcStateQ14[lpcHistoryIndex],
                decision.autoregressiveStateQ14[0],
                warpingQ16
            )
            val decisionAutoregressiveStateQ14 = decision.autoregressiveStateQ14
            var previousStateQ14 = smlawb(
                decisionAutoregressiveStateQ14[0],
                sub32(decisionAutoregressiveStateQ14[1], currentStateQ14),
                warpingQ16
            )
            decisionAutoregressiveStateQ14[0] = currentStateQ14
            var autoregressiveNoiseQ10 = smulwb(
                currentStateQ14,
                autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset]
            )
            var autoregressiveCoefficientIndex = 2
            while (autoregressiveCoefficientIndex < shapingLpcOrder) {
                currentStateQ14 = smlawb(
                    decisionAutoregressiveStateQ14[autoregressiveCoefficientIndex - 1],
                    sub32(
                        decisionAutoregressiveStateQ14[autoregressiveCoefficientIndex],
                        previousStateQ14
                    ),
                    warpingQ16
                )
                decisionAutoregressiveStateQ14[autoregressiveCoefficientIndex - 1] = previousStateQ14
                autoregressiveNoiseQ10 = smlawb(
                    autoregressiveNoiseQ10,
                    previousStateQ14,
                    autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset + autoregressiveCoefficientIndex - 1]
                )
                previousStateQ14 = smlawb(
                    decisionAutoregressiveStateQ14[autoregressiveCoefficientIndex],
                    sub32(
                        decisionAutoregressiveStateQ14[(autoregressiveCoefficientIndex + 1)],
                        currentStateQ14
                    ),
                    warpingQ16
                )
                decisionAutoregressiveStateQ14[autoregressiveCoefficientIndex] = currentStateQ14
                autoregressiveNoiseQ10 = smlawb(
                    autoregressiveNoiseQ10,
                    currentStateQ14,
                    autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset + autoregressiveCoefficientIndex]
                )
                autoregressiveCoefficientIndex += 2
            }
            decisionAutoregressiveStateQ14[shapingLpcOrder - 1] = previousStateQ14
            autoregressiveNoiseQ10 = smlawb(
                autoregressiveNoiseQ10,
                previousStateQ14,
                autoregressiveShapingCoefficientsQ13[autoregressiveShapingCoefficientOffset + shapingLpcOrder - 1]
            )
            autoregressiveNoiseQ10 = rightShift(autoregressiveNoiseQ10, 1)
            autoregressiveNoiseQ10 = smlawb(autoregressiveNoiseQ10, decision.lowFrequencyShapingStateQ12, tiltQ14)
            var lowFrequencyNoiseQ10 = leftShift(
                smulwb(decision.ltpShapeQ10[sampleHistoryIndex], lowFrequencyShapingQ14),
                2
            )
            lowFrequencyNoiseQ10 =
                smlawt(lowFrequencyNoiseQ10, decision.lowFrequencyShapingStateQ12, lowFrequencyShapingQ14)
            val shapingPredictionQ10 =
                sub32(
                    add32(rightShift(sub32(ltpPredQ14, longTermPredictionNoiseQ14), 4), lpcPredQ10),
                    add32(autoregressiveNoiseQ10, lowFrequencyNoiseQ10)
                )
            var residualQ10 = sub32(scaledInputQ10[sampleIndex], shapingPredictionQ10)
            residualQ10 = (residualQ10 xor dither) - dither
            residualQ10 = sub32(residualQ10, offsetQ10)
            residualQ10 = limit32(residualQ10, (-((64 shl 10))), (64 shl 10))
            var firstPulseQ10: Int
            var secondPulseQ10: Int
            var firstRateDistortionQ10: Int
            var secondRateDistortionQ10: Int
            if (residualQ10 < -1536) {
                firstPulseQ10 = leftShift(rshiftRound(residualQ10, 10), 10)
                residualQ10 = sub32(residualQ10, firstPulseQ10)
                firstRateDistortionQ10 = rightShift(
                    smlabb(
                        multiply((-add32(firstPulseQ10, offsetQ10)), lambdaQ10),
                        residualQ10,
                        residualQ10
                    ), 10
                )
                secondRateDistortionQ10 =
                    sub32(add32(firstRateDistortionQ10, 1024), addLShift32(lambdaQ10, residualQ10, 1))
                secondPulseQ10 = add32(firstPulseQ10, 1024)
            } else if (residualQ10 > 512) {
                firstPulseQ10 = leftShift(rshiftRound(residualQ10, 10), 10)
                residualQ10 = sub32(residualQ10, firstPulseQ10)
                firstRateDistortionQ10 = rightShift(
                    smlabb(
                        multiply(add32(firstPulseQ10, offsetQ10), lambdaQ10),
                        residualQ10,
                        residualQ10
                    ),
                    10
                )
                secondRateDistortionQ10 =
                    sub32(add32(firstRateDistortionQ10, 1024), subLShift32(lambdaQ10, residualQ10, 1))
                secondPulseQ10 = sub32(firstPulseQ10, 1024)
            } else {
                val offsetLambdaProductQ20 = smulbb(offsetQ10, lambdaQ10)
                secondRateDistortionQ10 = rightShift(
                    smlabb(offsetLambdaProductQ20, residualQ10, residualQ10),
                    10
                )
                firstRateDistortionQ10 = add32(secondRateDistortionQ10, 1024)
                firstRateDistortionQ10 = add32(
                    firstRateDistortionQ10,
                    subRShift32(addLShift32(lambdaQ10, residualQ10, 1), offsetLambdaProductQ20, 9)
                )
                firstPulseQ10 = -1024
                secondPulseQ10 = 0
            }
            if (firstRateDistortionQ10 < secondRateDistortionQ10) {
                bestState.rateDistortionQ10 = add32(decision.rateDistortionQ10, firstRateDistortionQ10)
                secondState.rateDistortionQ10 = add32(decision.rateDistortionQ10, secondRateDistortionQ10)
                bestState.pulseQ10 = firstPulseQ10
                secondState.pulseQ10 = secondPulseQ10
            } else {
                bestState.rateDistortionQ10 = add32(decision.rateDistortionQ10, secondRateDistortionQ10)
                secondState.rateDistortionQ10 = add32(decision.rateDistortionQ10, firstRateDistortionQ10)
                bestState.pulseQ10 = secondPulseQ10
                secondState.pulseQ10 = firstPulseQ10
            }
            var excitationQ10 = add32(offsetQ10, bestState.pulseQ10)
            excitationQ10 = (excitationQ10 xor dither) - dither
            var lpcExcitationQ10 = add32(excitationQ10, rshiftRound(ltpPredQ14, 4))
            var quantizedSignalQ10 = add32(lpcExcitationQ10, lpcPredQ10)
            var lowFrequencyShapingQ10 = sub32(quantizedSignalQ10, autoregressiveNoiseQ10)
            bestState.ltpShapeQ10 = sub32(lowFrequencyShapingQ10, lowFrequencyNoiseQ10)
            bestState.lowFrequencyShapingStateQ12 = leftShift(lowFrequencyShapingQ10, 2)
            bestState.quantizedSignalQ14 = leftShift(quantizedSignalQ10, 4)
            bestState.lpcExcitationQ16 = leftShift(lpcExcitationQ10, 6)
            excitationQ10 = add32(offsetQ10, secondState.pulseQ10)
            excitationQ10 = (excitationQ10 xor dither) - dither
            lpcExcitationQ10 = add32(excitationQ10, rshiftRound(ltpPredQ14, 4))
            quantizedSignalQ10 = add32(lpcExcitationQ10, lpcPredQ10)
            lowFrequencyShapingQ10 = sub32(quantizedSignalQ10, autoregressiveNoiseQ10)
            secondState.ltpShapeQ10 = sub32(lowFrequencyShapingQ10, lowFrequencyNoiseQ10)
            secondState.lowFrequencyShapingStateQ12 = leftShift(lowFrequencyShapingQ10, 2)
            secondState.quantizedSignalQ14 = leftShift(quantizedSignalQ10, 4)
            secondState.lpcExcitationQ16 = leftShift(lpcExcitationQ10, 6)
            decisionIndex++
        }
        sampleHistoryIndex = (sampleHistoryIndex - 1) and DECISION_DELAY_MASK
        val lastSampleIndex = (sampleHistoryIndex + decisionDelay) and DECISION_DELAY_MASK
        val winnerIndex = selectBestSampleStateIndex(sampleStates, delayedDecisionCount)
        penalizeMismatchedSeedStates(
            delayedDecisions,
            sampleStates,
            delayedDecisionCount,
            winnerIndex,
            lastSampleIndex
        )
        replaceWeakDelayedDecision(
            delayedDecisions,
            sampleStates,
            delayedDecisionCount,
            sampleIndex
        )
        val winner = delayedDecisions[winnerIndex]
        if (subframeIndex > 0 || sampleIndex >= decisionDelay) {
            quantizedPulses[((pulseOffset + sampleIndex) - decisionDelay)] =
                rightShift(winner.pulseQ10[lastSampleIndex], 10)
            quantizedSignalOutputQ10[((quantizedSignalOffset + sampleIndex) - decisionDelay)] =
                sat16(
                    rshiftRound(
                        smulww(winner.quantizedSignalQ10[lastSampleIndex], winner.gainQ16[lastSampleIndex]),
                        10
                    )
                )
            stateLtpShapeQ10[(ltpShapeWriteIndex - decisionDelay)] = winner.ltpShapeQ10[lastSampleIndex]
            scaledLtpExcitationQ16[(ltpExcitationWriteIndex - decisionDelay)] =
                winner.lpcExcitationQ16[lastSampleIndex]
        }
        ltpShapeWriteIndex++
        ltpExcitationWriteIndex++
        var historyUpdateIndex = 0
        while (historyUpdateIndex < delayedDecisionCount) {
            commitDelayedDecisionSample(
                delayedDecisions[historyUpdateIndex],
                sampleStates[historyUpdateIndex][0],
                sampleIndex,
                sampleHistoryIndex,
                gainQ16
            )
            historyUpdateIndex++
        }
        sampleIndex++
    }
    sampleHistoryIndexRef.value = sampleHistoryIndex
    state.ltpShapeWriteIndex = ltpShapeWriteIndex
    state.ltpExcitationWriteIndex = ltpExcitationWriteIndex
    resetDelayedDecisionHistories(delayedDecisions, delayedDecisionCount, subframeLength)
}

@Suppress("DuplicatedCode")
internal fun quantizeNoiseShape(
    context: NoiseShapeQuantizerContext,
    features: FrameFeatures,
    prefilteredInput: IntArray,
    quantizedPulses: IntArray,
    lsfInterpFactorQ2: Int,
    predCoefQ12: IntArray,
    ltpCoefQ14: IntArray,
    ar2Q13: IntArray,
    harmShapeGainQ14: IntArray,
    tiltQ14: IntArray,
    lowFrequencyShapingQ14: IntArray,
    gainsQ16: IntArray,
    lambdaQ10: Int,
    ltpScaleQ14: Int,
    useLbrr: Boolean = false,
) {
    val state = selectNoiseShapeQuantizerState(context, useLbrr)
    val frameLength = context.frameLength
    val subframeLength = context.subframeLength
    val predictionLpcOrder = context.predictionLpcOrder
    val shapingLpcOrder = context.shapingLpcOrder
    val offsetQ10 = quantizationOffsetsQ10[features.signalType][features.quantOffsetType]
    val lsfInterpolationFlag = if (lsfInterpFactorQ2 == 1 shl 2) 0 else 1
    val rewhiteningMask = rewhiteningSubframeMask(lsfInterpolationFlag)
    val scaledLtpExcitationQ16 = context.scaledLtpExcitationQ16
    val ltpExcitation = context.ltpExcitation
    val scaledInputQ10 = context.scaledInputQ10
    val lpcFilterState = context.lpcFilterState
    state.randomSeed = features.seed
    var lag = state.previousLag
    state.ltpShapeWriteIndex = frameLength
    state.ltpExcitationWriteIndex = frameLength
    var inputOffset = 0
    var pulseOffset = 0
    var quantizedSignalOffset = frameLength
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val predictionCoefficientOffset =
            ((subframeIndex shr 1) or (1 - lsfInterpolationFlag)) * MAX_LPC_ORDER
        val ltpCoefficientOffset = subframeIndex * LTP_ORDER
        val shapingCoefficientOffset = subframeIndex * MAX_SHAPE_LPC_ORDER
        val harmonicShapeGainQ14 = harmShapeGainQ14[subframeIndex]
        val harmShapeFirPackedQ14 = packHarmonicShapeFirQ14(harmonicShapeGainQ14)
        state.shouldRewhiten = false
        if (features.signalType == SIG_TYPE_VOICED) {
            lag = features.pitchL[subframeIndex]
            if ((subframeIndex and rewhiteningMask) == 0) {
                prepareRewhitening(
                    state,
                    frameLength,
                    lag,
                    predictionLpcOrder,
                    predCoefQ12,
                    predictionCoefficientOffset,
                    lpcFilterState,
                    ltpExcitation,
                    quantizedSignalStartIndex = frameLength - lag - predictionLpcOrder -
                            (LTP_ORDER shr 1) + subframeIndex * subframeLength
                )
            }
        }
        scaleNoiseShapeQuantizerState(
            state,
            prefilteredInput,
            inputOffset,
            scaledInputQ10,
            subframeLength,
            ltpExcitation,
            scaledLtpExcitationQ16,
            subframeIndex,
            ltpScaleQ14,
            gainsQ16,
            features.pitchL
        )
        quantizeNoiseShapeSubframe(
            state,
            features.signalType,
            scaledInputQ10,
            quantizedPulses,
            pulseOffset,
            state.quantizedSignalQ10,
            quantizedSignalOffset,
            scaledLtpExcitationQ16,
            predCoefQ12,
            predictionCoefficientOffset,
            ltpCoefQ14,
            ltpCoefficientOffset,
            ar2Q13,
            shapingCoefficientOffset,
            lag,
            harmShapeFirPackedQ14,
            tiltQ14[subframeIndex],
            lowFrequencyShapingQ14[subframeIndex],
            gainsQ16[subframeIndex],
            lambdaQ10,
            offsetQ10,
            subframeLength,
            shapingLpcOrder,
            predictionLpcOrder
        )
        inputOffset += subframeLength
        pulseOffset += subframeLength
        quantizedSignalOffset += subframeLength
        subframeIndex++
    }
    state.previousLag = features.pitchL[(NB_SUBFR - 1)]
    copyFinalNoiseShapeQuantizerState(state, frameLength)
}

@Suppress("DuplicatedCode")
internal fun quantizeNoiseShapeDelayedDecision(
    context: NoiseShapeQuantizerContext,
    features: FrameFeatures,
    prefilteredInput: IntArray,
    quantizedPulses: IntArray,
    lsfInterpFactorQ2: Int,
    predCoefQ12: IntArray,
    ltpCoefQ14: IntArray,
    ar2Q13: IntArray,
    harmShapeGainQ14: IntArray,
    tiltQ14: IntArray,
    lowFrequencyShapingQ14: IntArray,
    gainsQ16: IntArray,
    lambdaQ10: Int,
    ltpScaleQ14: Int,
    useLbrr: Boolean = false,
) {
    val state = selectNoiseShapeQuantizerState(context, useLbrr)
    val frameLength = context.frameLength
    val subframeLength = context.subframeLength
    val predictionLpcOrder = context.predictionLpcOrder
    val shapingLpcOrder = context.shapingLpcOrder
    val delayedDecisionCount = min(context.delayedDecisionCount, MAX_DELAYED_DECISION_STATES)
    val offsetQ10 = quantizationOffsetsQ10[features.signalType][features.quantOffsetType]
    val lsfInterpolationFlag = if (lsfInterpFactorQ2 == 1 shl 2) 0 else 1
    val rewhiteningMask = rewhiteningSubframeMask(lsfInterpolationFlag)
    val scaledLtpExcitationQ16 = context.scaledLtpExcitationQ16
    val ltpExcitation = context.ltpExcitation
    val scaledInputQ10 = context.scaledInputQ10
    val lpcFilterState = context.lpcFilterState
    val delayedDecisions = context.delayedDecisions
    val sampleStates = context.sampleStates
    var lag = state.previousLag
    val effectiveDecisionDelay = calculateEffectiveDecisionDelay(features, lag, subframeLength)
    initializeDelayedDecisionStates(
        state,
        delayedDecisions,
        delayedDecisionCount,
        frameLength,
        features.seed
    )
    val sampleHistoryIndexRef = context.sampleHistoryIndex
    sampleHistoryIndexRef.value = 0
    var inputOffset = 0
    var pulseOffset = 0
    var quantizedSignalOffset = frameLength
    var subframesSinceRewhitening = 0
    state.ltpShapeWriteIndex = frameLength
    state.ltpExcitationWriteIndex = frameLength
    var subframeIndex = 0
    while (subframeIndex < NB_SUBFR) {
        val predictionCoefficientOffset =
            ((subframeIndex shr 1) or (1 - lsfInterpolationFlag)) * MAX_LPC_ORDER
        val ltpCoefficientOffset = subframeIndex * LTP_ORDER
        val shapingCoefficientOffset = subframeIndex * MAX_SHAPE_LPC_ORDER
        val harmonicShapeGainQ14 = harmShapeGainQ14[subframeIndex]
        val harmShapeFirPackedQ14 = packHarmonicShapeFirQ14(harmonicShapeGainQ14)
        state.shouldRewhiten = false
        if (features.signalType == SIG_TYPE_VOICED) {
            lag = features.pitchL[subframeIndex]
            if ((subframeIndex and rewhiteningMask) == 0) {
                if (subframeIndex == 2) {
                    val winnerIndex =
                        selectMinimumDistortionIndex(delayedDecisions, delayedDecisionCount)
                    penalizeRejectedDecisions(delayedDecisions, delayedDecisionCount, winnerIndex)
                    val winner = delayedDecisions[winnerIndex]
                    flushDelayedDecisionSamples(
                        state,
                        winner,
                        quantizedPulses,
                        pulseOffset,
                        quantizedSignalOffset,
                        sampleHistoryIndexRef.value,
                        effectiveDecisionDelay
                    )
                    subframesSinceRewhitening = 0
                }
                prepareRewhitening(
                    state,
                    frameLength,
                    lag,
                    predictionLpcOrder,
                    predCoefQ12,
                    predictionCoefficientOffset,
                    lpcFilterState,
                    ltpExcitation,
                    quantizedSignalStartIndex = frameLength - lag - predictionLpcOrder -
                            (LTP_ORDER shr 1) + subframeIndex * subframeLength
                )
            }
        }
        scaleDelayedDecisionStates(
            state,
            delayedDecisions,
            prefilteredInput,
            inputOffset,
            scaledInputQ10,
            subframeLength,
            ltpExcitation,
            scaledLtpExcitationQ16,
            subframeIndex,
            delayedDecisionCount,
            ltpScaleQ14,
            gainsQ16,
            features.pitchL
        )
        quantizeNoiseShapeSubframeDelayedDecision(
            state,
            delayedDecisions,
            sampleStates,
            features.signalType,
            scaledInputQ10,
            quantizedPulses,
            pulseOffset,
            state.quantizedSignalQ10,
            quantizedSignalOffset,
            scaledLtpExcitationQ16,
            predCoefQ12,
            predictionCoefficientOffset,
            ltpCoefQ14,
            ltpCoefficientOffset,
            ar2Q13,
            shapingCoefficientOffset,
            lag,
            harmShapeFirPackedQ14,
            tiltQ14[subframeIndex],
            lowFrequencyShapingQ14[subframeIndex],
            gainsQ16[subframeIndex],
            lambdaQ10,
            offsetQ10,
            subframeLength,
            subframesSinceRewhitening,
            shapingLpcOrder,
            predictionLpcOrder,
            context.warpingQ16,
            delayedDecisionCount,
            sampleHistoryIndexRef,
            effectiveDecisionDelay
        )
        subframesSinceRewhitening++
        inputOffset += subframeLength
        pulseOffset += subframeLength
        quantizedSignalOffset += subframeLength
        subframeIndex++
    }
    val winnerIndex = selectMinimumDistortionIndex(delayedDecisions, delayedDecisionCount)
    val winner = delayedDecisions[winnerIndex]
    features.seed = winner.frameStartSeed
    flushDelayedDecisionSamples(
        state,
        winner,
        quantizedPulses,
        frameLength,
        2 * frameLength,
        sampleHistoryIndexRef.value,
        effectiveDecisionDelay
    )
    flushDelayedDecisionExcitation(
        state,
        winner,
        scaledLtpExcitationQ16,
        sampleHistoryIndexRef.value,
        effectiveDecisionDelay
    )
    copyWinnerQuantizerState(state, winner, subframeLength)
    state.previousLag = features.pitchL[(NB_SUBFR - 1)]
    copyFinalNoiseShapeQuantizerState(state, frameLength)
}
