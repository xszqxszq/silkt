package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.fixedpoint.*
import xyz.xszq.silkt.internal.model.MAX_LPC_ORDER
import xyz.xszq.silkt.internal.util.RefInt

// Ported from tsilk.
// Source: solve_LDL_FIX.ts

private const val MAX_MATRIX_SIZE = MAX_LPC_ORDER
private const val CONDITION_FACTOR_Q31 = 21474

internal class CholeskyWorkspace(
    val lowerTriangleQ16: IntArray,
    val intermediateSolution: IntArray,
    val workQ0: IntArray,
    val diagonalQ0: IntArray,
    val inverseDiagonalQ36: IntArray,
    val inverseDiagonalQ48: IntArray,
) {
    constructor() : this(
        lowerTriangleQ16 = IntArray(MAX_MATRIX_SIZE * MAX_MATRIX_SIZE),
        intermediateSolution = IntArray(MAX_MATRIX_SIZE),
        workQ0 = IntArray(MAX_MATRIX_SIZE),
        diagonalQ0 = IntArray(MAX_MATRIX_SIZE),
        inverseDiagonalQ36 = IntArray(MAX_MATRIX_SIZE),
        inverseDiagonalQ48 = IntArray(MAX_MATRIX_SIZE),
    )
}

private fun ldlFactorize(
    matrix: IntArray,
    matrixOffset: Int,
    matrixOrder: Int,
    lowerTriangleQ16: IntArray,
    workQ0: IntArray,
    diagonalQ0: IntArray,
    inverseDiagonalQ36: IntArray,
    inverseDiagonalQ48: IntArray,
) {
    var needsRegularization = true
    var round = 0
    val minimumDiagonal = max32(
        smmul(
            matrix[matrixOffset] + matrix[matrixOffset + matrixOrder * matrixOrder - 1],
            CONDITION_FACTOR_Q31,
        ),
        1 shl 9
    )

    while (needsRegularization && round < matrixOrder) {
        needsRegularization = false
        var column = 0
        while (column < matrixOrder) {
            val rowOffset = column * matrixOrder
            var dotProduct = 0
            var row = 0
            while (row < column) {
                workQ0[row] = smulww(diagonalQ0[row], lowerTriangleQ16[rowOffset + row])
                dotProduct = smlaww(dotProduct, workQ0[row], lowerTriangleQ16[rowOffset + row])
                row++
            }

            dotProduct = sub32(matrix[matrixOffset + column * matrixOrder + column], dotProduct)
            if (dotProduct < minimumDiagonal) {
                val regularizer = sub32(smulbb(round + 1, minimumDiagonal), dotProduct)
                var diagonalIndex = 0
                while (diagonalIndex < matrixOrder) {
                    val diagonalOffset = matrixOffset + diagonalIndex * matrixOrder + diagonalIndex
                    matrix[diagonalOffset] = add32(matrix[diagonalOffset], regularizer)
                    diagonalIndex++
                }
                needsRegularization = true
                break
            }

            diagonalQ0[column] = dotProduct
            val oneDivDiagonalQ36 = inverse32VarQ(dotProduct, 36)
            val oneDivDiagonalQ40 = leftShift(oneDivDiagonalQ36, 4)
            val residual = sub32(1 shl 24, smulww(dotProduct, oneDivDiagonalQ40))
            val oneDivDiagonalQ48 = smulww(residual, oneDivDiagonalQ40)
            inverseDiagonalQ36[column] = oneDivDiagonalQ36
            inverseDiagonalQ48[column] = oneDivDiagonalQ48
            lowerTriangleQ16[column * matrixOrder + column] = 65536

            var destinationOffset = (column + 1) * matrixOrder
            var destinationRow = column + 1
            while (destinationRow < matrixOrder) {
                var weightedSum = 0
                var sourceColumn = 0
                while (sourceColumn < column) {
                    weightedSum = smlaww(
                        weightedSum,
                        workQ0[sourceColumn],
                        lowerTriangleQ16[destinationOffset + sourceColumn]
                    )
                    sourceColumn++
                }
                weightedSum = sub32(matrix[matrixOffset + rowOffset + destinationRow], weightedSum)
                lowerTriangleQ16[destinationRow * matrixOrder + column] = add32(
                    smmul(weightedSum, oneDivDiagonalQ48),
                    rightShift(smulww(weightedSum, oneDivDiagonalQ36), 4)
                )
                destinationOffset += matrixOrder
                destinationRow++
            }
            column++
        }
        round++
    }
}

private fun divideByDiagonalsQ16(
    values: IntArray,
    inverseDiagonalQ36: IntArray,
    inverseDiagonalQ48: IntArray,
    matrixOrder: Int
) {
    var index = 0
    while (index < matrixOrder) {
        val value = values[index]
        val q36Part = inverseDiagonalQ36[index]
        val q48Part = inverseDiagonalQ48[index]
        values[index] = add32(smmul(value, q48Part), rightShift(smulww(value, q36Part), 4))
        index++
    }
}

private fun solveLowerTriangleQ16(
    lowerTriangleQ16: IntArray,
    matrixOrder: Int,
    rightHandSide: IntArray,
    solutionQ16: IntArray
) {
    var row = 0
    while (row < matrixOrder) {
        val rowOffset = row * matrixOrder
        var weightedSum = 0
        var column = 0
        while (column < row) {
            weightedSum = smlaww(
                weightedSum,
                lowerTriangleQ16[rowOffset + column],
                solutionQ16[column]
            )
            column++
        }
        solutionQ16[row] = sub32(rightHandSide[row], weightedSum)
        row++
    }
}

private fun solveUpperTriangleQ16(
    lowerTriangleQ16: IntArray,
    matrixOrder: Int,
    rightHandSide: IntArray,
    solutionQ16: IntArray
) {
    var row = matrixOrder - 1
    while (row >= 0) {
        var weightedSum = 0
        var column = matrixOrder - 1
        while (column > row) {
            weightedSum = smlaww(
                weightedSum,
                lowerTriangleQ16[column * matrixOrder + row],
                solutionQ16[column]
            )
            column--
        }
        solutionQ16[row] = sub32(rightHandSide[row], weightedSum)
        row--
    }
}

internal fun solveLDL(
    workspace: CholeskyWorkspace,
    matrix: IntArray,
    matrixOffset: Int,
    matrixOrder: Int,
    rightHandSide: IntArray,
    solutionQ16: IntArray,
) {
    ldlFactorize(
        matrix,
        matrixOffset,
        matrixOrder,
        workspace.lowerTriangleQ16,
        workspace.workQ0,
        workspace.diagonalQ0,
        workspace.inverseDiagonalQ36,
        workspace.inverseDiagonalQ48,
    )
    solveLowerTriangleQ16(
        workspace.lowerTriangleQ16,
        matrixOrder,
        rightHandSide,
        workspace.intermediateSolution,
    )
    divideByDiagonalsQ16(
        workspace.intermediateSolution,
        workspace.inverseDiagonalQ36,
        workspace.inverseDiagonalQ48,
        matrixOrder,
    )
    solveUpperTriangleQ16(
        workspace.lowerTriangleQ16,
        matrixOrder,
        workspace.intermediateSolution,
        solutionQ16,
    )
}

internal fun regularizeCorrelations(
    correlations: IntArray,
    correlationsOffset: Int,
    noiseTotal: RefInt,
    noise: Int,
    matrixOrder: Int
) {
    var diagonalIndex = 0
    while (diagonalIndex < matrixOrder) {
        val diagonalOffset = correlationsOffset + diagonalIndex * matrixOrder + diagonalIndex
        correlations[diagonalOffset] = add32(correlations[diagonalOffset], noise)
        diagonalIndex++
    }
    noiseTotal.value += noise
}
