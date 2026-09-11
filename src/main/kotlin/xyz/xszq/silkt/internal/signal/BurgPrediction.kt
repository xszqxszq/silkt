package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.util.RefInt

// Ported from tsilk.
// Source: burg_modified.ts

private const val PREDICTOR_QA = 25
private const val HEADROOM_BITS = 2
private const val MIN_RIGHT_SHIFTS = -16
private const val MAX_RIGHT_SHIFTS = 32 - PREDICTOR_QA

internal fun burgModified(
    residualEnergy: RefInt,
    residualEnergyQ: RefInt,
    coefficientsQ16: IntArray,
    samples: IntArray,
    samplesOffset: Int,
    subframeLength: Int,
    subframeCount: Int,
    whiteNoiseFractionQ32: Int,
    order: Int
) {
    val firstRowCorrelations = IntArray(MAX_LPC_ORDER)
    val lastRowCorrelations = IntArray(MAX_LPC_ORDER)
    val predictorsQA = IntArray(MAX_LPC_ORDER)
    val forwardCorrelations = IntArray(MAX_LPC_ORDER + 1)
    val backwardCorrelations = IntArray(MAX_LPC_ORDER + 1)

    val zeroLagOutput = RefInt(0)
    val shiftOutput = RefInt(0)
    sumSquaresShift(
        zeroLagOutput,
        shiftOutput,
        samples,
        samplesOffset,
        subframeCount * subframeLength
    )

    var zeroLagCorrelation = zeroLagOutput.value
    var rightShifts = shiftOutput.value
    var leadingZeros: Int
    var extraRightShifts: Int

    if (rightShifts > MAX_RIGHT_SHIFTS) {
        zeroLagCorrelation = lshift32(zeroLagCorrelation, rightShifts - MAX_RIGHT_SHIFTS)
        rightShifts = MAX_RIGHT_SHIFTS
    } else {
        leadingZeros = clz32(zeroLagCorrelation) - 1
        extraRightShifts = HEADROOM_BITS - leadingZeros
        if (extraRightShifts > 0) {
            extraRightShifts = min(extraRightShifts, MAX_RIGHT_SHIFTS - rightShifts)
            zeroLagCorrelation = rshift32(zeroLagCorrelation, extraRightShifts)
        } else {
            extraRightShifts = max(extraRightShifts, MIN_RIGHT_SHIFTS - rightShifts)
            zeroLagCorrelation = lshift32(zeroLagCorrelation, -extraRightShifts)
        }
        rightShifts += extraRightShifts
    }

    var subframeIndex = 0
    while (subframeIndex < subframeCount) {
        val subframeOffset = samplesOffset + subframeIndex * subframeLength
        var lag = 1
        while (lag < order + 1) {
            if (rightShifts > 0) {
                val product = innerProduct16(
                    samples,
                    subframeOffset,
                    samples,
                    subframeOffset + lag,
                    subframeLength - lag
                )
                firstRowCorrelations[lag - 1] += (product.toLong() shr rightShifts).toInt()
            } else {
                firstRowCorrelations[lag - 1] += lshift32(
                    innerProductAligned(
                        samples,
                        subframeOffset,
                        samples,
                        subframeOffset + lag,
                        subframeLength - lag
                    ),
                    -rightShifts
                )
            }
            lag++
        }
        subframeIndex++
    }

    firstRowCorrelations.copyInto(lastRowCorrelations)
    val regularizedZeroLag =
        zeroLagCorrelation + smmul(whiteNoiseFractionQ32, zeroLagCorrelation) + 1
    forwardCorrelations[0] = regularizedZeroLag
    backwardCorrelations[0] = regularizedZeroLag

    var orderIndex = 0
    var energy: Int
    var coefficientIndex: Int
    while (orderIndex < order) {
        if (rightShifts > -2) {
            var subframeIndex = 0
            while (subframeIndex < subframeCount) {
                val subframeOffset = samplesOffset + subframeIndex * subframeLength
                val leadingSample =
                    -lshift32(samples[subframeOffset + orderIndex], 16 - rightShifts)
                val trailingSample = -lshift32(
                    samples[subframeOffset + subframeLength - orderIndex - 1],
                    16 - rightShifts
                )
                var forwardPrediction =
                    lshift32(samples[subframeOffset + orderIndex], PREDICTOR_QA - 16)
                var backwardPrediction = lshift32(
                    samples[subframeOffset + subframeLength - orderIndex - 1],
                    PREDICTOR_QA - 16
                )

                coefficientIndex = 0
                while (coefficientIndex < orderIndex) {
                    val predictorQA = predictorsQA[coefficientIndex]
                    val leadingHistoryIndex = subframeOffset + orderIndex - coefficientIndex - 1
                    val trailingHistoryIndex =
                        subframeOffset + subframeLength - orderIndex + coefficientIndex

                    firstRowCorrelations[coefficientIndex] = smlawb(
                        firstRowCorrelations[coefficientIndex],
                        leadingSample,
                        samples[leadingHistoryIndex]
                    )
                    lastRowCorrelations[coefficientIndex] = smlawb(
                        lastRowCorrelations[coefficientIndex],
                        trailingSample,
                        samples[trailingHistoryIndex]
                    )
                    forwardPrediction = smlawb(
                        forwardPrediction,
                        predictorQA,
                        samples[leadingHistoryIndex]
                    )
                    backwardPrediction = smlawb(
                        backwardPrediction,
                        predictorQA,
                        samples[trailingHistoryIndex]
                    )
                    coefficientIndex++
                }

                forwardPrediction = lshift32(
                    -forwardPrediction,
                    32 - PREDICTOR_QA - rightShifts
                )
                backwardPrediction = lshift32(
                    -backwardPrediction,
                    32 - PREDICTOR_QA - rightShifts
                )

                coefficientIndex = 0
                while (coefficientIndex <= orderIndex) {
                    forwardCorrelations[coefficientIndex] = smlawb(
                        forwardCorrelations[coefficientIndex],
                        forwardPrediction,
                        samples[subframeOffset + orderIndex - coefficientIndex]
                    )
                    backwardCorrelations[coefficientIndex] = smlawb(
                        backwardCorrelations[coefficientIndex],
                        backwardPrediction,
                        samples[subframeOffset + subframeLength - orderIndex + coefficientIndex - 1]
                    )
                    coefficientIndex++
                }
                subframeIndex++
            }
        } else {
            var subframeIndex = 0
            while (subframeIndex < subframeCount) {
                val subframeOffset = samplesOffset + subframeIndex * subframeLength
                val leadingSample =
                    -lshift32(samples[subframeOffset + orderIndex], -rightShifts)
                val trailingSample = -lshift32(
                    samples[subframeOffset + subframeLength - orderIndex - 1],
                    -rightShifts
                )
                var forwardPrediction = lshift32(samples[subframeOffset + orderIndex], 17)
                var backwardPrediction = lshift32(
                    samples[subframeOffset + subframeLength - orderIndex - 1],
                    17
                )

                coefficientIndex = 0
                while (coefficientIndex < orderIndex) {
                    val scaledPredictor =
                        rshiftRound(predictorsQA[coefficientIndex], PREDICTOR_QA - 17)
                    val leadingHistoryIndex = subframeOffset + orderIndex - coefficientIndex - 1
                    val trailingHistoryIndex =
                        subframeOffset + subframeLength - orderIndex + coefficientIndex

                    firstRowCorrelations[coefficientIndex] = mla(
                        firstRowCorrelations[coefficientIndex],
                        leadingSample,
                        samples[leadingHistoryIndex]
                    )
                    lastRowCorrelations[coefficientIndex] = mla(
                        lastRowCorrelations[coefficientIndex],
                        trailingSample,
                        samples[trailingHistoryIndex]
                    )
                    forwardPrediction = mla(
                        forwardPrediction,
                        samples[leadingHistoryIndex],
                        scaledPredictor
                    )
                    backwardPrediction = mla(
                        backwardPrediction,
                        samples[trailingHistoryIndex],
                        scaledPredictor
                    )
                    coefficientIndex++
                }

                forwardPrediction = -forwardPrediction
                backwardPrediction = -backwardPrediction

                coefficientIndex = 0
                while (coefficientIndex <= orderIndex) {
                    forwardCorrelations[coefficientIndex] = smlaww(
                        forwardCorrelations[coefficientIndex],
                        forwardPrediction,
                        lshift32(
                            samples[subframeOffset + orderIndex - coefficientIndex],
                            -rightShifts - 1
                        )
                    )
                    backwardCorrelations[coefficientIndex] = smlaww(
                        backwardCorrelations[coefficientIndex],
                        backwardPrediction,
                        lshift32(
                            samples[
                                subframeOffset + subframeLength - orderIndex + coefficientIndex - 1
                            ],
                            -rightShifts - 1
                        )
                    )
                    coefficientIndex++
                }
                subframeIndex++
            }
        }

        var forwardCorrelation = firstRowCorrelations[orderIndex]
        var backwardCorrelation = lastRowCorrelations[orderIndex]
        var numerator = 0
        energy = backwardCorrelations[0] + forwardCorrelations[0]

        coefficientIndex = 0
        while (coefficientIndex < orderIndex) {
            val predictorQA = predictorsQA[coefficientIndex]
            leadingZeros = min(32 - PREDICTOR_QA, clz32(abs(predictorQA)) - 1)
            val shiftedPredictor = lshift32(predictorQA, leadingZeros)
            val correlationIndex = orderIndex - coefficientIndex - 1
            val outputShift = 32 - PREDICTOR_QA - leadingZeros

            forwardCorrelation = addLshift32(
                forwardCorrelation,
                smmul(lastRowCorrelations[correlationIndex], shiftedPredictor),
                outputShift
            )
            backwardCorrelation = addLshift32(
                backwardCorrelation,
                smmul(firstRowCorrelations[correlationIndex], shiftedPredictor),
                outputShift
            )
            numerator = addLshift32(
                numerator,
                smmul(backwardCorrelations[orderIndex - coefficientIndex], shiftedPredictor),
                outputShift
            )
            energy = addLshift32(
                energy,
                smmul(
                    backwardCorrelations[coefficientIndex + 1] +
                            forwardCorrelations[coefficientIndex + 1],
                    shiftedPredictor
                ),
                outputShift
            )
            coefficientIndex++
        }

        forwardCorrelations[orderIndex + 1] = forwardCorrelation
        backwardCorrelations[orderIndex + 1] = backwardCorrelation
        numerator += backwardCorrelation
        numerator = lshift32(-numerator, 1)

        if (abs(numerator) >= energy) {
            var remainingIndex = orderIndex
            while (remainingIndex < order) {
                predictorsQA[remainingIndex] = 0
                remainingIndex++
            }
            break
        }

        val reflectionQ31 = div32VarQ(numerator, energy, 31)

        coefficientIndex = 0
        while (coefficientIndex < (orderIndex + 1 shr 1)) {
            val firstPredictor = predictorsQA[coefficientIndex]
            val secondPredictor = predictorsQA[orderIndex - coefficientIndex - 1]
            predictorsQA[coefficientIndex] = addLshift32(
                firstPredictor,
                smmul(secondPredictor, reflectionQ31),
                1
            )
            predictorsQA[orderIndex - coefficientIndex - 1] = addLshift32(
                secondPredictor,
                smmul(firstPredictor, reflectionQ31),
                1
            )
            coefficientIndex++
        }
        predictorsQA[orderIndex] = rshift32(reflectionQ31, 31 - PREDICTOR_QA)

        coefficientIndex = 0
        while (coefficientIndex <= orderIndex + 1) {
            val forwardValue = forwardCorrelations[coefficientIndex]
            val backwardValue = backwardCorrelations[orderIndex - coefficientIndex + 1]
            forwardCorrelations[coefficientIndex] = addLshift32(
                forwardValue,
                smmul(backwardValue, reflectionQ31),
                1
            )
            backwardCorrelations[orderIndex - coefficientIndex + 1] = addLshift32(
                backwardValue,
                smmul(forwardValue, reflectionQ31),
                1
            )
            coefficientIndex++
        }
        orderIndex++
    }

    energy = forwardCorrelations[0]
    var gainQ16 = 1 shl 16
    coefficientIndex = 0
    while (coefficientIndex < order) {
        val coefficientQ16 = rshiftRound(predictorsQA[coefficientIndex], PREDICTOR_QA - 16)
        energy = smlaww(energy, forwardCorrelations[coefficientIndex + 1], coefficientQ16)
        gainQ16 = smlaww(gainQ16, coefficientQ16, coefficientQ16)
        coefficientsQ16[coefficientIndex] = -coefficientQ16
        coefficientIndex++
    }

    residualEnergy.value = smlaww(
        energy,
        smmul(whiteNoiseFractionQ32, zeroLagCorrelation),
        -gainQ16
    )
    residualEnergyQ.value = -rightShifts
}
