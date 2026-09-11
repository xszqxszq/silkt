package xyz.xszq.silkt.internal.model

internal const val SIG_TYPE_VOICED: Int = 0
internal const val SIG_TYPE_UNVOICED: Int = 1
internal const val FRAME_LENGTH_MS: Int = 20
internal const val MAX_FS_KHZ: Int = 24
internal const val MAX_FRAME_LENGTH: Int = FRAME_LENGTH_MS * MAX_FS_KHZ
internal const val MAX_LPC_ORDER: Int = 16
internal const val MAX_SHAPE_LPC_ORDER: Int = 16
internal const val MAX_SHAPE_LPC_WIN_LENGTH: Int = 480
internal const val LTP_ORDER: Int = 5
internal const val NB_SUBFR: Int = 4
internal const val MAX_ARITHM_BYTES: Int = 1024
