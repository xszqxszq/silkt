package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.NB_SUBFR
import xyz.xszq.silkt.internal.quantization.NlsfConversionWorkspace
import xyz.xszq.silkt.internal.quantization.convertNlsfToStablePrediction
import xyz.xszq.silkt.internal.quantization.interpolateValues
import xyz.xszq.silkt.internal.quantization.predictionCoefficientsToNlsf
import xyz.xszq.silkt.internal.util.RefInt


private const val TUNING_FIND_LPC_COND_FAC: Double = 0.000025
private const val TUNING_FIND_LPC_CHIRP: Double = 0.99995

// Ported from tsilk.
// Source: corrMatrix_FIX.ts

internal fun corrVector(
    input: IntArray,
    inputOffset: Int,
    target: IntArray,
    targetOffset: Int,
    windowLength: Int,
    order: Int,
    correlations: IntArray,
    rightShifts: Int
) {
    var inputWindowEndOffset = inputOffset + order - 1

    if (rightShifts > 0) {
        var lag = 0
        while (lag < order) {
            var innerProduct = 0
            var sampleIndex = 0
            while (sampleIndex < windowLength) {
                innerProduct += rshift32(
                    smulbb(input[inputWindowEndOffset + sampleIndex], target[targetOffset + sampleIndex]),
                    rightShifts
                )
                sampleIndex++
            }
            correlations[lag] = innerProduct
            inputWindowEndOffset--
            lag++
        }
    } else {
        var lag = 0
        while (lag < order) {
            correlations[lag] = innerProductAligned(input, inputWindowEndOffset, target, targetOffset, windowLength)
            inputWindowEndOffset--
            lag++
        }
    }
}

internal fun corrMatrix(
    energyOutput: RefInt,
    shiftOutput: RefInt,
    input: IntArray,
    inputOffset: Int,
    windowLength: Int,
    order: Int,
    headRoom: Int,
    correlationMatrix: IntArray,
    matrixOffset: Int,
    outputRightShifts: RefInt
) {
    sumSquaresShift(energyOutput, shiftOutput, input, inputOffset, windowLength + order - 1)

    var energy = energyOutput.value
    var rightShifts = shiftOutput.value
    val headRoomRightShifts = max(headRoom - clz32(energy), 0)
    energy = rshift32(energy, headRoomRightShifts)
    rightShifts += headRoomRightShifts

    var sampleIndex = 0
    while (sampleIndex < order - 1) {
        energy -= rshift32(smulbb(input[inputOffset + sampleIndex], input[inputOffset + sampleIndex]), rightShifts)
        sampleIndex++
    }

    if (rightShifts < outputRightShifts.value) {
        energy = rshift32(energy, outputRightShifts.value - rightShifts)
        rightShifts = outputRightShifts.value
    }

    correlationMatrix[matrixOffset] = energy
    val firstWindowEndOffset = inputOffset + order - 1
    var diagonalIndex = 1
    while (diagonalIndex < order) {
        val futureEndOffset = firstWindowEndOffset + windowLength - diagonalIndex
        val pastEndOffset = firstWindowEndOffset - diagonalIndex
        energy = sub32(energy, rshift32(smulbb(input[futureEndOffset], input[futureEndOffset]), rightShifts))
        energy = add32(energy, rshift32(smulbb(input[pastEndOffset], input[pastEndOffset]), rightShifts))
        correlationMatrix[matrixOffset + diagonalIndex * order + diagonalIndex] = energy
        diagonalIndex++
    }

    var laggedWindowEndOffset = inputOffset + order - 2
    if (rightShifts > 0) {
        var lag = 1
        while (lag < order) {
            energy = 0
            var sampleOffset = 0
            while (sampleOffset < windowLength) {
                energy += rshift32(
                    smulbb(
                        input[firstWindowEndOffset + sampleOffset],
                        input[laggedWindowEndOffset + sampleOffset]
                    ),
                    rightShifts
                )
                sampleOffset++
            }

            correlationMatrix[matrixOffset + lag * order] = energy
            correlationMatrix[matrixOffset + lag] = energy

            var lagOffset = 1
            while (lagOffset < order - lag) {
                val firstFutureOffset = firstWindowEndOffset + windowLength - lagOffset
                val laggedFutureOffset = laggedWindowEndOffset + windowLength - lagOffset
                val firstPastOffset = firstWindowEndOffset - lagOffset
                val laggedPastOffset = laggedWindowEndOffset - lagOffset
                energy = sub32(
                    energy,
                    rshift32(smulbb(input[firstFutureOffset], input[laggedFutureOffset]), rightShifts)
                )
                energy = add32(
                    energy,
                    rshift32(smulbb(input[firstPastOffset], input[laggedPastOffset]), rightShifts)
                )
                correlationMatrix[matrixOffset + (lag + lagOffset) * order + lagOffset] = energy
                correlationMatrix[matrixOffset + lagOffset * order + lag + lagOffset] = energy
                lagOffset++
            }
            laggedWindowEndOffset--
            lag++
        }
    } else {
        var lag = 1
        while (lag < order) {
            energy = innerProductAligned(input, firstWindowEndOffset, input, laggedWindowEndOffset, windowLength)
            correlationMatrix[matrixOffset + lag * order] = energy
            correlationMatrix[matrixOffset + lag] = energy

            var lagOffset = 1
            while (lagOffset < order - lag) {
                val firstFutureOffset = firstWindowEndOffset + windowLength - lagOffset
                val laggedFutureOffset = laggedWindowEndOffset + windowLength - lagOffset
                val firstPastOffset = firstWindowEndOffset - lagOffset
                val laggedPastOffset = laggedWindowEndOffset - lagOffset
                energy = sub32(energy, smulbb(input[firstFutureOffset], input[laggedFutureOffset]))
                energy = smlabb(energy, input[firstPastOffset], input[laggedPastOffset])
                correlationMatrix[matrixOffset + (lag + lagOffset) * order + lagOffset] = energy
                correlationMatrix[matrixOffset + lagOffset * order + lag + lagOffset] = energy
                lagOffset++
            }
            laggedWindowEndOffset--
            lag++
        }
    }
    outputRightShifts.value = rightShifts
}

// Ported from tsilk.
// Source: find_Lpc_FIX.ts

private val findLpcConditionFactor: Int = fixConst(TUNING_FIND_LPC_COND_FAC, 32)

private val findLpcChirp: Int = fixConst(TUNING_FIND_LPC_CHIRP, 16)

@Suppress("DuplicatedCode")
internal fun findLpc(
    nlsfQ15: IntArray,
    interpIndex: RefInt,
    prevNlsfQ15: IntArray,
    workspace: LinearPredictionWorkspace,
    conversionWorkspace: NlsfConversionWorkspace,
    useInterpolatedNlsfs: Int,
    lpcOrder: Int,
    input: IntArray,
    inputOffset: Int,
    subfrLength: Int
) {
    val coefficientsQ16 = workspace.coefficientsQ16
    val tailCoefficientsQ16 = workspace.tailCoefficientsQ16
    val candidateCoefficientsQ12 = workspace.candidateCoefficientsQ12
    val interpolatedNlsfQ15 = workspace.interpolatedNlsfQ15
    val filterState = workspace.filterState
    val residual = workspace.residual

    val firstHalfEnergy = workspace.firstHalfEnergyOutput
    val secondHalfEnergy = workspace.secondHalfEnergyOutput
    val firstHalfShift = workspace.firstHalfShiftOutput
    val secondHalfShift = workspace.secondHalfShiftOutput
    val residualEnergyOutput = workspace.residualEnergyOutput
    val residualEnergyQOutput = workspace.residualEnergyQOutput
    val tailEnergyOutput = workspace.tailEnergyOutput
    val tailEnergyQOutput = workspace.tailEnergyQOutput

    var shift: Int
    var interpolatedEnergy: Int
    var interpolatedEnergyQ: Int
    var residualEnergy: Int
    var residualEnergyQ: Int
    var tailEnergy: Int
    var tailEnergyQ: Int

    interpIndex.value = 4
    burgModified(
        residualEnergyOutput,
        residualEnergyQOutput,
        coefficientsQ16,
        input,
        inputOffset,
        subfrLength,
        NB_SUBFR,
        findLpcConditionFactor,
        lpcOrder
    )
    residualEnergy = residualEnergyOutput.value
    residualEnergyQ = residualEnergyQOutput.value
    expandBandwidth32(coefficientsQ16, lpcOrder, findLpcChirp)

    if (useInterpolatedNlsfs == 1) {
        burgModified(
            tailEnergyOutput,
            tailEnergyQOutput,
            tailCoefficientsQ16,
            input,
            inputOffset + (NB_SUBFR shr 1) * subfrLength,
            subfrLength,
            NB_SUBFR shr 1,
            findLpcConditionFactor,
            lpcOrder
        )
        tailEnergy = tailEnergyOutput.value
        tailEnergyQ = tailEnergyQOutput.value
        expandBandwidth32(tailCoefficientsQ16, lpcOrder, findLpcChirp)

        shift = tailEnergyQ - residualEnergyQ
        if (shift >= 0) {
            if ((shift < 32)) {
                residualEnergy -= rightShift(tailEnergy, shift)
            }
        } else {
            residualEnergy = rightShift(residualEnergy, -shift) - tailEnergy
            residualEnergyQ = tailEnergyQ
        }

        predictionCoefficientsToNlsf(nlsfQ15, tailCoefficientsQ16, lpcOrder, conversionWorkspace)
        var interpolationIndex = 3
        while (interpolationIndex >= 0) {
            interpolateValues(
                interpolatedNlsfQ15,
                prevNlsfQ15,
                nlsfQ15,
                interpolationIndex,
                lpcOrder
            )
            convertNlsfToStablePrediction(
                candidateCoefficientsQ12,
                interpolatedNlsfQ15,
                lpcOrder,
                conversionWorkspace
            )
            filterState.fill(0, 0, lpcOrder)
            lpcAnalysisFilter(
                input,
                inputOffset,
                candidateCoefficientsQ12,
                filterState,
                residual,
                2 * subfrLength,
                lpcOrder
            )
            sumSquaresShift(
                firstHalfEnergy,
                firstHalfShift,
                residual,
                lpcOrder,
                subfrLength - lpcOrder
            )
            sumSquaresShift(
                secondHalfEnergy,
                secondHalfShift,
                residual,
                lpcOrder + subfrLength,
                subfrLength - lpcOrder,
            )
            var firstEnergy = firstHalfEnergy.value
            var secondEnergy = secondHalfEnergy.value
            shift = firstHalfShift.value - secondHalfShift.value
            if (shift >= 0) {
                secondEnergy = rightShift(secondEnergy, shift)
                interpolatedEnergyQ = -firstHalfShift.value
            } else {
                firstEnergy = rightShift(firstEnergy, -shift)
                interpolatedEnergyQ = -secondHalfShift.value
            }
            interpolatedEnergy = add32(firstEnergy, secondEnergy)

            shift = interpolatedEnergyQ - residualEnergyQ
            val interpolatedIsLower = if (shift >= 0) {
                rightShift(interpolatedEnergy, shift) < residualEnergy
            } else if (-shift < 32) {
                interpolatedEnergy < rightShift(residualEnergy, -shift)
            } else {
                false
            }

            if (interpolatedIsLower) {
                residualEnergy = interpolatedEnergy
                residualEnergyQ = interpolatedEnergyQ
                interpIndex.value = interpolationIndex
            }
            interpolationIndex--
        }
    }

    if (interpIndex.value == 4) {
        predictionCoefficientsToNlsf(nlsfQ15, coefficientsQ16, lpcOrder, conversionWorkspace)
    }
}
