package xyz.xszq.silkt.internal.tables

// Ported from tsilk.
// Source: tables/tables_other.ts

internal object SharedEncoderTables {

    internal val samplingRatesTable: IntArray = intArrayOf(8, 12, 16, 24)

    internal val samplingRatesCdf: IntArray = intArrayOf(0, 16000, 32000, 48000, 65535)

    internal val typeOffsetCdf: IntArray = intArrayOf(0, 37522, 41030, 44212, 65535)

    internal val typeOffsetJointCdf: Array<IntArray> = arrayOf(
        intArrayOf(0, 57686, 61230, 62358, 65535), intArrayOf(0, 18346, 40067, 43659, 65535),
        intArrayOf(0, 22694, 24279, 35507, 65535), intArrayOf(0, 6067, 7215, 13010, 65535),
    )

    internal val nlsfInterpolationFactorCdf: IntArray = intArrayOf(0, 3706, 8703, 19226, 30926, 65535)

    internal val frameTerminationCdf: IntArray = intArrayOf(0, 20000, 45000, 56000, 65535)

    internal val seedCdf: IntArray = intArrayOf(0, 16384, 32768, 49152, 65535)

    internal val voiceActivityCdf: IntArray = intArrayOf(0, 22000, 65535)

    internal val lsbCdf: IntArray = intArrayOf(0, 40000, 65535)

    internal val ltpScaleCdf: IntArray = intArrayOf(0, 32000, 48000, 65535)

    internal val ltpScalesTableQ14: IntArray = intArrayOf(15565, 11469, 8192)

    internal val targetRateTableNb: IntArray = intArrayOf(0, 8000, 9000, 11000, 13000, 16000, 22000, 100000)

    internal val targetRateTableMb: IntArray = intArrayOf(0, 10000, 12000, 14000, 17000, 21000, 28000, 100000)

    internal val targetRateTableWb: IntArray = intArrayOf(0, 11000, 14000, 17000, 21000, 26000, 36000, 100000)

    internal val targetRateTableSwb: IntArray = intArrayOf(0, 13000, 16000, 19000, 25000, 32000, 46000, 100000)

    internal val snrTableQ1: IntArray = intArrayOf(19, 31, 35, 39, 43, 47, 54, 64)

    internal val signCdf: IntArray = intArrayOf(
        37840, 36944, 36251, 35304, 34715, 35503, 34529, 34296, 34016, 47659, 44945, 42503, 40235, 38569, 40254, 37851,
        37243, 36595, 43410, 44121, 43127, 40978, 38845, 40433, 38252, 37795, 36637, 59159, 55630, 51806, 48073,
        45036, 48416, 43857, 42678, 41146,
    )

    internal val transitionLpBQ28: Array<IntArray> = arrayOf(
        intArrayOf(250767114, 501534038, 250767114), intArrayOf(209867381, 419732057, 209867381),
        intArrayOf(170987846, 341967853, 170987846), intArrayOf(131531482, 263046905, 131531482),
        intArrayOf(89306658, 178584282, 89306658),
    )

    internal val transitionLpAQ28: Array<IntArray> = arrayOf(
        intArrayOf(506393414, 239854379), intArrayOf(411067935, 169683996), intArrayOf(306733530, 116694253),
        intArrayOf(185807084, 77959395), intArrayOf(35497197, 57401098),
    )

}
