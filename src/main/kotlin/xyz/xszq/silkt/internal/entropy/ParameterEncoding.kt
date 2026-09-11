package xyz.xszq.silkt.internal.entropy

import xyz.xszq.silkt.internal.fixedpoint.abs
import xyz.xszq.silkt.internal.fixedpoint.rightShift
import xyz.xszq.silkt.internal.model.MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.model.SIG_TYPE_VOICED
import xyz.xszq.silkt.internal.quantization.NlsfCodebook
import xyz.xszq.silkt.internal.signal.FrameFeatures
import xyz.xszq.silkt.internal.tables.ExcitationTables.deltaGainCdf
import xyz.xszq.silkt.internal.tables.ExcitationTables.gainCdf
import xyz.xszq.silkt.internal.tables.ExcitationTables.maxPulsesTable
import xyz.xszq.silkt.internal.tables.ExcitationTables.pulsesPerBlockBitsQ6
import xyz.xszq.silkt.internal.tables.ExcitationTables.pulsesPerBlockCdf
import xyz.xszq.silkt.internal.tables.ExcitationTables.rateLevelsBitsQ6
import xyz.xszq.silkt.internal.tables.ExcitationTables.rateLevelsCdf
import xyz.xszq.silkt.internal.tables.ExcitationTables.shellCodeTable0
import xyz.xszq.silkt.internal.tables.ExcitationTables.shellCodeTable1
import xyz.xszq.silkt.internal.tables.ExcitationTables.shellCodeTable2
import xyz.xszq.silkt.internal.tables.ExcitationTables.shellCodeTable3
import xyz.xszq.silkt.internal.tables.ExcitationTables.shellCodeTableOffsets
import xyz.xszq.silkt.internal.tables.LtpTables.ltpGainCdfTables
import xyz.xszq.silkt.internal.tables.LtpTables.ltpPerIndexCdf
import xyz.xszq.silkt.internal.tables.PitchTables.pitchContourCdf
import xyz.xszq.silkt.internal.tables.PitchTables.pitchContourNbCdf
import xyz.xszq.silkt.internal.tables.PitchTables.pitchLagMbCdf
import xyz.xszq.silkt.internal.tables.PitchTables.pitchLagNbCdf
import xyz.xszq.silkt.internal.tables.PitchTables.pitchLagSwbCdf
import xyz.xszq.silkt.internal.tables.PitchTables.pitchLagWbCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.lsbCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.ltpScaleCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.nlsfInterpolationFactorCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.samplingRatesCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.samplingRatesTable
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.seedCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.signCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.typeOffsetCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.typeOffsetJointCdf
import xyz.xszq.silkt.internal.tables.SharedEncoderTables.voiceActivityCdf
import xyz.xszq.silkt.internal.util.RefInt

private const val SHELL_CODEC_FRAME_LENGTH = 16
private const val MAX_NB_SHELL_BLOCKS = MAX_FRAME_LENGTH / SHELL_CODEC_FRAME_LENGTH
private const val N_RATE_LEVELS = 10
private const val MAX_PULSES = 18

internal class PulseEncodingWorkspace {
    val absolutePulses = IntArray(MAX_FRAME_LENGTH)
    val blockPulseSums = IntArray(MAX_NB_SHELL_BLOCKS)
    val rightShiftCounts = IntArray(MAX_NB_SHELL_BLOCKS)
    val combinedPulses = IntArray(8)
    val pairSums = IntArray(8)
    val quadSums = IntArray(4)
    val twoSums = IntArray(2)
    val totalSums = IntArray(1)
    val signCdf = IntArray(3)
}

// Ported from tsilk.
// Source: encode_pulses.ts

@Suppress("SameParameterValue")
private fun combineAndCheck(
    combinedPulses: IntArray,
    inputPulses: IntArray,
    maxPulses: Int,
    pairCount: Int,
    inputOffset: Int,
    outputOffset: Int
): Int {
    var pairIndex = 0
    while (pairIndex < pairCount) {
        val sum = inputPulses[inputOffset + 2 * pairIndex] +
                inputPulses[inputOffset + 2 * pairIndex + 1]
        if (sum > maxPulses) {
            return 1
        }
        combinedPulses[outputOffset + pairIndex] = sum
        pairIndex++
    }
    return 0
}

private fun encodeSplit(
    rangeCoder: RangeCoderState,
    childPulseCount: Int,
    totalPulseCount: Int,
    shellCdfTable: IntArray
) {
    if (totalPulseCount > 0) {
        rangeEncode(
            rangeCoder,
            childPulseCount,
            shellCdfTable,
            shellCodeTableOffsets[totalPulseCount]
        )
    }
}

private fun shellEncode(
    rangeCoder: RangeCoderState,
    workspace: PulseEncodingWorkspace,
    blockPulses: IntArray,
    blockOffset: Int
) {
    combineAndCheck(workspace.pairSums, blockPulses, Int.MAX_VALUE, 8, blockOffset, 0)
    combineAndCheck(workspace.quadSums, workspace.pairSums, Int.MAX_VALUE, 4, 0, 0)
    combineAndCheck(workspace.twoSums, workspace.quadSums, Int.MAX_VALUE, 2, 0, 0)
    combineAndCheck(workspace.totalSums, workspace.twoSums, Int.MAX_VALUE, 1, 0, 0)

    encodeSplit(rangeCoder, workspace.twoSums[0], workspace.totalSums[0], shellCodeTable3)
    encodeSplit(rangeCoder, workspace.quadSums[0], workspace.twoSums[0], shellCodeTable2)
    encodeSplit(rangeCoder, workspace.pairSums[0], workspace.quadSums[0], shellCodeTable1)
    encodeSplit(rangeCoder, blockPulses[blockOffset], workspace.pairSums[0], shellCodeTable0)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 2], workspace.pairSums[1], shellCodeTable0)
    encodeSplit(rangeCoder, workspace.pairSums[2], workspace.quadSums[1], shellCodeTable1)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 4], workspace.pairSums[2], shellCodeTable0)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 6], workspace.pairSums[3], shellCodeTable0)
    encodeSplit(rangeCoder, workspace.quadSums[2], workspace.twoSums[1], shellCodeTable2)
    encodeSplit(rangeCoder, workspace.pairSums[4], workspace.quadSums[2], shellCodeTable1)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 8], workspace.pairSums[4], shellCodeTable0)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 10], workspace.pairSums[5], shellCodeTable0)
    encodeSplit(rangeCoder, workspace.pairSums[6], workspace.quadSums[3], shellCodeTable1)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 12], workspace.pairSums[6], shellCodeTable0)
    encodeSplit(rangeCoder, blockPulses[blockOffset + 14], workspace.pairSums[7], shellCodeTable0)
}

private fun encodeSigns(
    rangeCoder: RangeCoderState,
    pulses: IntArray,
    length: Int,
    signalType: Int,
    quantOffsetType: Int,
    rateLevelIndex: Int,
    workspace: PulseEncodingWorkspace
) {
    val signCdfIndex =
        (N_RATE_LEVELS - 1) * ((signalType shl 1) + quantOffsetType) + rateLevelIndex
    workspace.signCdf[1] = signCdf[signCdfIndex]
    workspace.signCdf[2] = 65535
    var pulseIndex = 0
    while (pulseIndex < length) {
        if (pulses[pulseIndex] != 0) {
            val negativeSign = rightShift(pulses[pulseIndex], 15) + 1
            rangeEncode(rangeCoder, negativeSign, workspace.signCdf)
        }
        pulseIndex++
    }
}

private fun encodePulses(
    rangeCoder: RangeCoderState,
    signalType: Int,
    quantOffsetType: Int,
    pulses: IntArray,
    frameLength: Int,
    workspace: PulseEncodingWorkspace
) {
    val shellBlockCount = frameLength / SHELL_CODEC_FRAME_LENGTH
    normalizePulseBlocks(pulses, frameLength, shellBlockCount, workspace)
    val rateLevelIndex = selectRateLevelIndex(signalType, shellBlockCount, workspace)

    rangeEncode(rangeCoder, rateLevelIndex, rateLevelsCdf[signalType])
    encodePulseBlockSums(rangeCoder, shellBlockCount, rateLevelIndex, workspace)
    encodeShellBlocks(rangeCoder, shellBlockCount, workspace)
    encodePulseShiftBits(rangeCoder, pulses, shellBlockCount, workspace)

    encodeSigns(
        rangeCoder,
        pulses,
        frameLength,
        signalType,
        quantOffsetType,
        rateLevelIndex,
        workspace
    )
}

private fun normalizePulseBlocks(
    pulses: IntArray,
    frameLength: Int,
    shellBlockCount: Int,
    workspace: PulseEncodingWorkspace
) {
    val absolutePulses = workspace.absolutePulses
    val blockPulseSums = workspace.blockPulseSums
    val rightShiftCounts = workspace.rightShiftCounts
    val combinedPulses = workspace.combinedPulses

    var sampleIndex = 0
    while (sampleIndex < frameLength) {
        absolutePulses[sampleIndex] = abs(pulses[sampleIndex])
        sampleIndex++
    }

    var blockOffset = 0
    var shellBlockIndex = 0
    while (shellBlockIndex < shellBlockCount) {
        rightShiftCounts[shellBlockIndex] = 0
        while (true) {
            var overflowCount = 0
            overflowCount += combineAndCheck(
                combinedPulses,
                absolutePulses,
                maxPulsesTable[0],
                8,
                blockOffset,
                0
            )
            overflowCount += combineAndCheck(
                combinedPulses,
                combinedPulses,
                maxPulsesTable[1],
                4,
                0,
                0
            )
            overflowCount += combineAndCheck(
                combinedPulses,
                combinedPulses,
                maxPulsesTable[2],
                2,
                0,
                0
            )

            blockPulseSums[shellBlockIndex] = combinedPulses[0] + combinedPulses[1]
            if (blockPulseSums[shellBlockIndex] > maxPulsesTable[3]) {
                overflowCount++
            }

            if (overflowCount == 0) {
                break
            }

            rightShiftCounts[shellBlockIndex]++
            var sampleInBlockIndex = 0
            while (sampleInBlockIndex < SHELL_CODEC_FRAME_LENGTH) {
                val sampleOffset = blockOffset + sampleInBlockIndex
                absolutePulses[sampleOffset] = rightShift(absolutePulses[sampleOffset], 1)
                sampleInBlockIndex++
            }
        }
        blockOffset += SHELL_CODEC_FRAME_LENGTH
        shellBlockIndex++
    }
}

private fun selectRateLevelIndex(
    signalType: Int,
    shellBlockCount: Int,
    workspace: PulseEncodingWorkspace
): Int {
    val blockPulseSums = workspace.blockPulseSums
    val rightShiftCounts = workspace.rightShiftCounts
    var rateLevelIndex = 0
    var minimumSumBitsQ6 = Int.MAX_VALUE
    var rateLevelCandidate = 0
    while (rateLevelCandidate < N_RATE_LEVELS - 1) {
        val blockBitsQ6 = pulsesPerBlockBitsQ6[rateLevelCandidate]
        var sumBitsQ6 = rateLevelsBitsQ6[signalType][rateLevelCandidate]
        var blockIndex = 0
        while (blockIndex < shellBlockCount) {
            sumBitsQ6 += if (rightShiftCounts[blockIndex] > 0) {
                blockBitsQ6[MAX_PULSES + 1]
            } else {
                blockBitsQ6[blockPulseSums[blockIndex]]
            }
            blockIndex++
        }
        if (sumBitsQ6 < minimumSumBitsQ6) {
            minimumSumBitsQ6 = sumBitsQ6
            rateLevelIndex = rateLevelCandidate
        }
        rateLevelCandidate++
    }
    return rateLevelIndex
}

private fun encodePulseBlockSums(
    rangeCoder: RangeCoderState,
    shellBlockCount: Int,
    rateLevelIndex: Int,
    workspace: PulseEncodingWorkspace
) {
    val blockPulseSums = workspace.blockPulseSums
    val rightShiftCounts = workspace.rightShiftCounts
    val blockCdf = pulsesPerBlockCdf[rateLevelIndex]
    var encodedBlockIndex = 0
    while (encodedBlockIndex < shellBlockCount) {
        if (rightShiftCounts[encodedBlockIndex] == 0) {
            rangeEncode(rangeCoder, blockPulseSums[encodedBlockIndex], blockCdf)
        } else {
            rangeEncode(rangeCoder, MAX_PULSES + 1, blockCdf)
            var escapedBlockCount = 0
            while (escapedBlockCount < rightShiftCounts[encodedBlockIndex] - 1) {
                rangeEncode(rangeCoder, MAX_PULSES + 1, pulsesPerBlockCdf[N_RATE_LEVELS - 1])
                escapedBlockCount++
            }
            rangeEncode(
                rangeCoder,
                blockPulseSums[encodedBlockIndex],
                pulsesPerBlockCdf[N_RATE_LEVELS - 1]
            )
        }
        encodedBlockIndex++
    }
}

private fun encodeShellBlocks(
    rangeCoder: RangeCoderState,
    shellBlockCount: Int,
    workspace: PulseEncodingWorkspace
) {
    val absolutePulses = workspace.absolutePulses
    val blockPulseSums = workspace.blockPulseSums
    var shellBlockIndex = 0
    while (shellBlockIndex < shellBlockCount) {
        if (blockPulseSums[shellBlockIndex] > 0) {
            shellEncode(
                rangeCoder,
                workspace,
                absolutePulses,
                shellBlockIndex * SHELL_CODEC_FRAME_LENGTH
            )
        }
        shellBlockIndex++
    }
}

private fun encodePulseShiftBits(
    rangeCoder: RangeCoderState,
    pulses: IntArray,
    shellBlockCount: Int,
    workspace: PulseEncodingWorkspace
) {
    val rightShiftCounts = workspace.rightShiftCounts
    var lsbBlockIndex = 0
    while (lsbBlockIndex < shellBlockCount) {
        if (rightShiftCounts[lsbBlockIndex] > 0) {
            val lsbBlockOffset = lsbBlockIndex * SHELL_CODEC_FRAME_LENGTH
            val remainingShiftCount = rightShiftCounts[lsbBlockIndex] - 1
            var sampleInBlockIndex = 0
            while (sampleInBlockIndex < SHELL_CODEC_FRAME_LENGTH) {
                val sampleOffset = lsbBlockOffset + sampleInBlockIndex
                // The source truncates the absolute pulse to signed 8 bits.
                val wrappedAbsolutePulse = (abs(pulses[sampleOffset]) shl 24) shr 24
                var bitIndex = remainingShiftCount
                while (bitIndex > 0) {
                    rangeEncode(
                        rangeCoder,
                        (wrappedAbsolutePulse shr bitIndex) and 1,
                        lsbCdf
                    )
                    bitIndex--
                }
                rangeEncode(rangeCoder, wrappedAbsolutePulse and 1, lsbCdf)
                sampleInBlockIndex++
            }
        }
        lsbBlockIndex++
    }
}

// Ported from tsilk.
// Source: encode_parameters.ts

private fun encodeNlsfIndices(
    rangeCoder: RangeCoderState,
    codebook: NlsfCodebook,
    nlsfIndices: IntArray
) {
    var stageIndex = 0
    while (stageIndex < codebook.stageCount) {
        val cdfStartOffset = codebook.cdfStartOffsets[stageIndex]
        rangeEncode(
            rangeCoder,
            nlsfIndices[stageIndex],
            codebook.cdf,
            cdfStartOffset
        )
        stageIndex++
    }
}

internal fun encodeParameters(
    payloadFrameCount: Int,
    samplingRateKHz: Int,
    typeOffsetPrevious: RefInt,
    nlsfCodebooks: Array<NlsfCodebook>,
    frameLength: Int,
    voiceActivityFlag: Boolean,
    features: FrameFeatures,
    rangeCoder: RangeCoderState,
    pulses: IntArray,
    pulseWorkspace: PulseEncodingWorkspace
) {
    if (payloadFrameCount == 0) {
        var samplingRateIndex = 0
        while (samplingRateIndex < samplingRatesTable.size) {
            if (samplingRatesTable[samplingRateIndex] == samplingRateKHz) {
                break
            }
            samplingRateIndex++
        }
        rangeEncode(rangeCoder, samplingRateIndex, samplingRatesCdf)
    }

    val typeOffset = 2 * features.signalType + features.quantOffsetType
    if (payloadFrameCount == 0) {
        rangeEncode(rangeCoder, typeOffset, typeOffsetCdf)
    } else {
        rangeEncode(rangeCoder, typeOffset, typeOffsetJointCdf[typeOffsetPrevious.value])
    }
    typeOffsetPrevious.value = typeOffset

    if (payloadFrameCount == 0) {
        rangeEncode(rangeCoder, features.gainsIndices[0], gainCdf[features.signalType])
    } else {
        rangeEncode(rangeCoder, features.gainsIndices[0], deltaGainCdf)
    }
    var subframeIndex = 1
    while (subframeIndex < NB_SUBFR) {
        rangeEncode(rangeCoder, features.gainsIndices[subframeIndex], deltaGainCdf)
        subframeIndex++
    }

    val nlsfCodebook = nlsfCodebooks[features.signalType]
    encodeNlsfIndices(rangeCoder, nlsfCodebook, features.nlsfIndices)
    rangeEncode(rangeCoder, features.nlsfInterpCoefQ2, nlsfInterpolationFactorCdf)

    if (features.signalType == SIG_TYPE_VOICED) {
        val pitchLagCdf = when (samplingRateKHz) {
            8 -> pitchLagNbCdf
            12 -> pitchLagMbCdf
            16 -> pitchLagWbCdf
            else -> pitchLagSwbCdf
        }
        rangeEncode(rangeCoder, features.lagIndex, pitchLagCdf)

        val pitchContourCdf =
            if (samplingRateKHz == 8) pitchContourNbCdf else pitchContourCdf
        rangeEncode(rangeCoder, features.contourIndex, pitchContourCdf)

        rangeEncode(rangeCoder, features.perIndex, ltpPerIndexCdf)
        var ltpSubframeIndex = 0
        while (ltpSubframeIndex < NB_SUBFR) {
            rangeEncode(
                rangeCoder,
                features.ltpIndex[ltpSubframeIndex],
                ltpGainCdfTables[features.perIndex]
            )
            ltpSubframeIndex++
        }
        rangeEncode(rangeCoder, features.ltpScaleIndex, ltpScaleCdf)
    }

    rangeEncode(rangeCoder, features.seed, seedCdf)
    encodePulses(
        rangeCoder,
        features.signalType,
        features.quantOffsetType,
        pulses,
        frameLength,
        pulseWorkspace
    )
    rangeEncode(rangeCoder, if (voiceActivityFlag) 1 else 0, voiceActivityCdf)
}
