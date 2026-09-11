package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.entropy.RangeCoderState
import xyz.xszq.silkt.internal.model.MAX_ARITHM_BYTES

internal class EncoderLbrrState {
    var oldestPacketIndex: Int = 0
    var previousLastGainIndex: Int = 0
    val rangeCoder: RangeCoderState = RangeCoderState()
    val packets: Array<LbrrPacket> =
        Array(MAX_LBRR_DELAY) { LbrrPacket(payload = IntArray(MAX_ARITHM_BYTES)) }

    fun resetPackets() {
        for (packet in packets) {
            packet.byteCount = 0
            packet.usage = SILK_NO_LBRR
        }
    }

    fun resetHistory() {
        oldestPacketIndex = 0
        resetPackets()
    }
}
