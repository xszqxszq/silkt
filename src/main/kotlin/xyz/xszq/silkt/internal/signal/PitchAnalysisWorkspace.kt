package xyz.xszq.silkt.internal.signal

import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_D_SRCH_LENGTH
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MAX_DECIMATE_STATE_LENGTH
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MAX_FRAME_LENGTH
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MAX_FRAME_LENGTH_ST1
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MAX_FRAME_LENGTH_ST2
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_MAX_LAG
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_CBKS_STAGE_2_EXT
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_CBKS_STAGE_3_MAX
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_STAGE_3_LAGS
import xyz.xszq.silkt.internal.tables.PitchTables.PITCH_EST_NB_SUB_FR
import xyz.xszq.silkt.internal.util.RefInt

private const val STAGE3_LAG_VALUE_COUNT = 22
private const val FIND_PITCH_LPC_WIN_MAX = 576

/** Reusable pitch-analysis storage; each stage overwrites the parts it reads. */
internal class PitchAnalysisWorkspace {
    val signal8KHz = IntArray(PITCH_EST_MAX_FRAME_LENGTH_ST2)
    val signal4KHz = IntArray(PITCH_EST_MAX_FRAME_LENGTH_ST1)
    val analysisSignal = IntArray(PITCH_EST_MAX_FRAME_LENGTH)
    val decimationState = IntArray(PITCH_EST_MAX_DECIMATE_STATE_LENGTH)
    val correlations = Array(PITCH_EST_NB_SUB_FR) { IntArray((PITCH_EST_MAX_LAG shr 1) + 5) }
    val searchLags = IntArray(PITCH_EST_D_SRCH_LENGTH)
    val lagCandidates = IntArray((PITCH_EST_MAX_LAG shr 1) + 5)
    val stage2Scores = IntArray(PITCH_EST_NB_CBKS_STAGE_2_EXT)
    val stage3Correlations =
        IntArray(PITCH_EST_NB_SUB_FR * PITCH_EST_NB_CBKS_STAGE_3_MAX * PITCH_EST_NB_STAGE_3_LAGS)
    val stage3Energies =
        IntArray(PITCH_EST_NB_SUB_FR * PITCH_EST_NB_CBKS_STAGE_3_MAX * PITCH_EST_NB_STAGE_3_LAGS)
    val stage3LagValues = IntArray(STAGE3_LAG_VALUE_COUNT)
    val weightedSignal = IntArray(FIND_PITCH_LPC_WIN_MAX)
    val autoCorrelation = IntArray(MAX_FIND_PITCH_LPC_ORDER + 1)
    val reflectionCoefficientsQ15 = IntArray(MAX_FIND_PITCH_LPC_ORDER)
    val coefficientsQ24 = IntArray(MAX_FIND_PITCH_LPC_ORDER)
    val coefficientsQ12 = IntArray(MAX_FIND_PITCH_LPC_ORDER)
    val predictionState = IntArray(MAX_FIND_PITCH_LPC_ORDER)
    val autocorrelationScale = RefInt()
    val lagAndContourIndices = IntArray(2)
    val ltpCorrelationOutput = IntArray(1)
}
