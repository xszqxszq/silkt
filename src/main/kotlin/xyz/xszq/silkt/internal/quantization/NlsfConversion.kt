package xyz.xszq.silkt.internal.quantization

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.signal.expandBandwidth
import xyz.xszq.silkt.internal.signal.expandBandwidth32
import xyz.xszq.silkt.internal.signal.lpcInversePredGain
import xyz.xszq.silkt.internal.tables.NlsfTables.LSF_COSINE_TABLE_SIZE
import xyz.xszq.silkt.internal.tables.NlsfTables.lsfCosineTableQ12

private const val MAX_LPC_STABILIZE_ITERATIONS = 20
private const val BINARY_SEARCH_STEPS = 3
private const val POLYNOMIAL_Q16 = 16
private const val MAX_ROOT_SEARCH_RESTARTS = 30
private const val MAX_NLSF_PREDICTION_ORDER = 16
private const val NLSF_COSINE_INDEX_SHIFT = 8

internal class NlsfConversionWorkspace(
    val polynomialP: IntArray,
    val polynomialQ: IntArray,
    val cosineNlsfQ20: IntArray,
    val coefficientsQ12: IntArray,
) {
    val polynomials: Array<IntArray> = arrayOf(polynomialP, polynomialQ)

    constructor() : this(
        polynomialP = IntArray((MAX_LPC_ORDER shr 1) + 1),
        polynomialQ = IntArray((MAX_LPC_ORDER shr 1) + 1),
        cosineNlsfQ20 = IntArray(MAX_NLSF_PREDICTION_ORDER),
        coefficientsQ12 = IntArray(MAX_NLSF_PREDICTION_ORDER),
    )
}

// Ported from tsilk.
// Source: A2Nlsf.ts

private fun transformPolynomial(polynomial: IntArray, halfOrder: Int) {
    var exponent = 2
    while (exponent <= halfOrder) {
        var coefficientIndex = halfOrder
        while (coefficientIndex > exponent) {
            polynomial[coefficientIndex - 2] -= polynomial[coefficientIndex]
            coefficientIndex--
        }
        polynomial[exponent - 2] -= leftShift(polynomial[exponent], 1)
        exponent++
    }
}

private fun evaluatePolynomial(
    polynomial: IntArray,
    cosineQ12: Int,
    halfOrder: Int
): Int {
    var valueQ16 = polynomial[halfOrder]
    val cosineQ16 = leftShift(cosineQ12, 4)

    var coefficientIndex = halfOrder - 1
    while (coefficientIndex >= 0) {
        valueQ16 = smlaww(polynomial[coefficientIndex], valueQ16, cosineQ16)
        coefficientIndex--
    }
    return valueQ16
}

private fun initializeNlsfPolynomials(
    predictionCoefficientsQ16: IntArray,
    polynomialP: IntArray,
    polynomialQ: IntArray,
    halfOrder: Int
) {
    polynomialP[halfOrder] = leftShift(1, POLYNOMIAL_Q16)
    polynomialQ[halfOrder] = leftShift(1, POLYNOMIAL_Q16)

    var coefficientIndex = 0
    while (coefficientIndex < halfOrder) {
        polynomialP[coefficientIndex] =
            -predictionCoefficientsQ16[halfOrder - coefficientIndex - 1] -
                    predictionCoefficientsQ16[halfOrder + coefficientIndex]
        polynomialQ[coefficientIndex] =
            -predictionCoefficientsQ16[halfOrder - coefficientIndex - 1] +
                    predictionCoefficientsQ16[halfOrder + coefficientIndex]
        coefficientIndex++
    }

    coefficientIndex = halfOrder
    while (coefficientIndex > 0) {
        polynomialP[coefficientIndex - 1] -= polynomialP[coefficientIndex]
        polynomialQ[coefficientIndex - 1] += polynomialQ[coefficientIndex]
        coefficientIndex--
    }

    transformPolynomial(polynomialP, halfOrder)
    transformPolynomial(polynomialQ, halfOrder)
}

internal fun predictionCoefficientsToNlsf(
    nlsfQ15: IntArray,
    predictionCoefficientsQ16: IntArray,
    predictionOrder: Int,
    workspace: NlsfConversionWorkspace
) {
    val halfOrder = predictionOrder shr 1
    val polynomialP = workspace.polynomialP
    val polynomialQ = workspace.polynomialQ
    val polynomialPQ = workspace.polynomials

    initializeNlsfPolynomials(
        predictionCoefficientsQ16,
        polynomialP,
        polynomialQ,
        halfOrder
    )

    var polynomial = polynomialP
    var lowerCosine = lsfCosineTableQ12[0]
    var lowerValue = evaluatePolynomial(polynomial, lowerCosine, halfOrder)
    var rootIndex = 0
    if (lowerValue < 0) {
        nlsfQ15[0] = 0
        polynomial = polynomialQ
        lowerValue = evaluatePolynomial(polynomial, lowerCosine, halfOrder)
        rootIndex = 1
    }

    var cosineIndex = 1
    var restartCount = 0
    while (true) {
        var upperCosine = lsfCosineTableQ12[cosineIndex]
        var upperValue = evaluatePolynomial(polynomial, upperCosine, halfOrder)
        if (
            (lowerValue <= 0 && upperValue >= 0) ||
            (lowerValue >= 0 && upperValue <= 0)
        ) {
            var rootFraction = -256
            var searchStep = 0
            while (searchStep < BINARY_SEARCH_STEPS) {
                val middleCosine = rshiftRound(lowerCosine + upperCosine, 1)
                val middleValue = evaluatePolynomial(polynomial, middleCosine, halfOrder)
                if (
                    (lowerValue <= 0 && middleValue >= 0) ||
                    (lowerValue >= 0 && middleValue <= 0)
                ) {
                    upperCosine = middleCosine
                    upperValue = middleValue
                } else {
                    lowerCosine = middleCosine
                    lowerValue = middleValue
                    rootFraction = addRshift(rootFraction, 128, searchStep)
                }
                searchStep++
            }

            if (abs(lowerValue) < 65536) {
                val denominator = lowerValue - upperValue
                val numerator =
                    leftShift(lowerValue, 8 - BINARY_SEARCH_STEPS) +
                            rightShift(denominator, 1)
                if (denominator != 0) {
                    rootFraction += div32(numerator, denominator)
                }
            } else {
                rootFraction += div32(
                    lowerValue,
                    rightShift(lowerValue - upperValue, 8 - BINARY_SEARCH_STEPS)
                )
            }

            nlsfQ15[rootIndex] =
                min(lshift32(cosineIndex, 8) + rootFraction, 32767)
            rootIndex++
            if (rootIndex >= predictionOrder) {
                break
            }

            polynomial = polynomialPQ[rootIndex and 1]
            lowerCosine = lsfCosineTableQ12[cosineIndex - 1]
            lowerValue = leftShift(1 - (rootIndex and 2), 12)
        } else {
            cosineIndex++
            lowerCosine = upperCosine
            lowerValue = upperValue
            if (cosineIndex > LSF_COSINE_TABLE_SIZE) {
                restartCount++
                if (restartCount > MAX_ROOT_SEARCH_RESTARTS) {
                    fillUniformNlsf(nlsfQ15, predictionOrder)
                    return
                }

                expandBandwidth32(
                    predictionCoefficientsQ16,
                    predictionOrder,
                    65536 - smulbb(10 + restartCount, restartCount)
                )
                initializeNlsfPolynomials(
                    predictionCoefficientsQ16,
                    polynomialP,
                    polynomialQ,
                    halfOrder
                )

                polynomial = polynomialP
                lowerCosine = lsfCosineTableQ12[0]
                lowerValue = evaluatePolynomial(polynomial, lowerCosine, halfOrder)
                rootIndex = 0
                if (lowerValue < 0) {
                    nlsfQ15[0] = 0
                    polynomial = polynomialQ
                    lowerValue = evaluatePolynomial(polynomial, lowerCosine, halfOrder)
                    rootIndex = 1
                }
                cosineIndex = 1
            }
        }
    }
}

// Ported from tsilk.
// Source: nlsf2a_stable.ts

internal fun convertNlsfToStablePrediction(
    predictionCoefficientsQ12: IntArray,
    nlsfQ15: IntArray,
    predictionOrder: Int,
    workspace: NlsfConversionWorkspace
): Int {
    convertNlsfToPrediction(predictionCoefficientsQ12, nlsfQ15, predictionOrder, workspace)

    var stableIterations = 0
    var iteration = 0
    while (iteration < MAX_LPC_STABILIZE_ITERATIONS) {
        val unstable = lpcInversePredGain(
            predictionCoefficientsQ12,
            predictionOrder
        ).unstable
        if (!unstable) {
            stableIterations = iteration
            break
        }

        val chirpQ16 = 65536 - smulbb(10 + iteration, iteration)
        expandBandwidth(predictionCoefficientsQ12, 0, predictionOrder, chirpQ16)
        stableIterations = iteration + 1
        iteration++
    }
    return stableIterations
}

// Ported from tsilk.
// Source: nlsf_a.ts

@Suppress("SameParameterValue")
private fun findNlsfPolynomial(
    output: IntArray,
    outputOffset: Int,
    cosineNlsfQ20: IntArray,
    cosineOffset: Int,
    halfOrder: Int
) {
    output[outputOffset] = 1 shl 20
    output[outputOffset + 1] = -cosineNlsfQ20[cosineOffset]

    var order = 1
    while (order < halfOrder) {
        val cosineQ20 = cosineNlsfQ20[cosineOffset + 2 * order]
        val leadingProduct = cosineQ20.toLong() * output[outputOffset + order].toLong()
        output[outputOffset + order + 1] = toInt32(
            (output[outputOffset + order - 1] shl 1) -
                    ((leadingProduct + (1L shl 19)) shr 20).toInt()
        )

        var coefficientIndex = order
        while (coefficientIndex > 1) {
            val product =
                cosineQ20.toLong() * output[outputOffset + coefficientIndex - 1].toLong()
            output[outputOffset + coefficientIndex] = toInt32(
                output[outputOffset + coefficientIndex] +
                        output[outputOffset + coefficientIndex - 2] -
                        ((product + (1L shl 19)) shr 20).toInt()
            )
            coefficientIndex--
        }
        output[outputOffset + 1] -= cosineQ20
        order++
    }
}

private fun convertNlsfToPrediction(
    predictionCoefficientsQ12: IntArray,
    nlsfQ15: IntArray,
    predictionOrder: Int,
    workspace: NlsfConversionWorkspace
) {
    val cosineNlsfQ20 = workspace.cosineNlsfQ20
    val polynomialP = workspace.polynomialP
    val polynomialQ = workspace.polynomialQ
    val coefficientsQ12 = workspace.coefficientsQ12

    initializeNlsfCosines(nlsfQ15, cosineNlsfQ20, predictionOrder)

    val halfOrder = predictionOrder shr 1
    findNlsfPolynomial(polynomialP, 0, cosineNlsfQ20, 0, halfOrder)
    findNlsfPolynomial(polynomialQ, 0, cosineNlsfQ20, 1, halfOrder)
    combineNlsfPolynomials(
        coefficientsQ12,
        polynomialP,
        polynomialQ,
        predictionOrder,
        halfOrder
    )
    limitNlsfPredictionMagnitude(coefficientsQ12, predictionOrder)
    writeNlsfPredictionCoefficients(
        predictionCoefficientsQ12,
        coefficientsQ12,
        predictionOrder
    )
}

private fun initializeNlsfCosines(
    nlsfQ15: IntArray,
    cosineNlsfQ20: IntArray,
    predictionOrder: Int
) {
    var coefficientIndex = 0
    while (coefficientIndex < predictionOrder) {
        val cosineIndex = nlsfQ15[coefficientIndex] shr NLSF_COSINE_INDEX_SHIFT
        val cosineFraction =
            nlsfQ15[coefficientIndex] - (cosineIndex shl NLSF_COSINE_INDEX_SHIFT)
        val cosineQ12 = lsfCosineTableQ12[cosineIndex]
        val cosineDeltaQ12 = lsfCosineTableQ12[cosineIndex + 1] - cosineQ12
        cosineNlsfQ20[coefficientIndex] = toInt32(
            (cosineQ12 shl 8) + multiply(cosineDeltaQ12, cosineFraction)
        )
        coefficientIndex++
    }
}

private fun combineNlsfPolynomials(
    coefficientsQ12: IntArray,
    polynomialP: IntArray,
    polynomialQ: IntArray,
    predictionOrder: Int,
    halfOrder: Int
) {
    var coefficientIndex = 0
    while (coefficientIndex < halfOrder) {
        val pValue = toInt32(polynomialP[coefficientIndex + 1] + polynomialP[coefficientIndex])
        val qValue = toInt32(polynomialQ[coefficientIndex + 1] - polynomialQ[coefficientIndex])
        coefficientsQ12[coefficientIndex] =
            -rshiftRound(toInt32(pValue + qValue), 9)
        coefficientsQ12[predictionOrder - coefficientIndex - 1] =
            rshiftRound(toInt32(qValue - pValue), 9)
        coefficientIndex++
    }
}

private fun limitNlsfPredictionMagnitude(
    coefficientsQ12: IntArray,
    predictionOrder: Int
) {
    var scalingIteration = 0
    while (scalingIteration < 10) {
        var maximumMagnitude = 0
        var maximumIndex = 0
        var coefficientIndex = 0

        while (coefficientIndex < predictionOrder) {
            val magnitude = abs(coefficientsQ12[coefficientIndex])
            if (magnitude > maximumMagnitude) {
                maximumMagnitude = magnitude
                maximumIndex = coefficientIndex
            }
            coefficientIndex++
        }

        if (maximumMagnitude <= 32767) {
            break
        }

        maximumMagnitude = min(maximumMagnitude, 98369)
        val scalingQ16 = toInt32(
            65470 - div32(
                multiply(65470 shr 2, maximumMagnitude - 32767),
                rightShift(multiply(maximumMagnitude, maximumIndex + 1), 2)
            )
        )
        applyNlsfPredictionChirp(coefficientsQ12, predictionOrder, scalingQ16)
        scalingIteration++
    }
}

private fun writeNlsfPredictionCoefficients(
    predictionCoefficientsQ12: IntArray,
    coefficientsQ12: IntArray,
    predictionOrder: Int
) {
    var coefficientIndex = 0
    while (coefficientIndex < predictionOrder) {
        predictionCoefficientsQ12[coefficientIndex] =
            toInt16(sat16(coefficientsQ12[coefficientIndex]))
        coefficientIndex++
    }
}

@Suppress("DuplicatedCode")
private fun applyNlsfPredictionChirp(
    coefficientsQ12: IntArray,
    predictionOrder: Int,
    chirpQ16: Int
) {
    var chirp = chirpQ16

    var coefficientIndex = 0
    while (coefficientIndex < predictionOrder - 1) {
        coefficientsQ12[coefficientIndex] =
            toInt32(smulww(chirp, coefficientsQ12[coefficientIndex]))
        chirp = toInt32(smulww(chirp, chirpQ16))
        coefficientIndex++
    }
    coefficientsQ12[predictionOrder - 1] =
        toInt32(smulww(chirp, coefficientsQ12[predictionOrder - 1]))
}

private fun fillUniformNlsf(nlsfQ15: IntArray, predictionOrder: Int) {
    nlsfQ15[0] = div3216(1 shl 15, predictionOrder + 1)
    var coefficientIndex = 1
    while (coefficientIndex < predictionOrder) {
        nlsfQ15[coefficientIndex] = smulbb(coefficientIndex + 1, nlsfQ15[0])
        coefficientIndex++
    }
}
