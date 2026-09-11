package xyz.xszq.silkt

import xyz.xszq.silkt.internal.codec.EncoderOptions
import xyz.xszq.silkt.internal.codec.EncoderState
import xyz.xszq.silkt.internal.codec.SilkEncoder
import xyz.xszq.silkt.internal.codec.encodeFrame
import xyz.xszq.silkt.internal.model.MAX_ARITHM_BYTES
import xyz.xszq.silkt.internal.signal.resample
import xyz.xszq.silkt.internal.util.RefInt
import java.io.ByteArrayOutputStream
import java.io.InputStream
import java.io.OutputStream

/** Pure Kotlin SILK v3 encoder with a silk-codec compatible container. */
object SilkCoder {
    private val silkHeader = byteArrayOf(0x23, 0x21, 0x53, 0x49, 0x4c, 0x4b, 0x5f, 0x56, 0x33)
    private val supportedSampleRates = intArrayOf(8000, 12000, 16000, 24000)

    fun encode(
        pcm: ByteArray,
        sampleRate: Int = 24000,
        bitRate: Int = 24000,
        tencent: Boolean = true,
    ): ByteArray {
        val silk = ByteArrayOutputStream(16 + pcm.size / 3)
        encode(pcm, silk, sampleRate, bitRate, tencent)
        return silk.toByteArray()
    }

    fun encode(
        pcm: InputStream,
        silk: OutputStream,
        sampleRate: Int = 24000,
        bitRate: Int = 24000,
    ) {
        encode(readAll(pcm), silk, sampleRate, bitRate, true)
    }

    fun encode(
        pcm: InputStream,
        silk: OutputStream,
        tencent: Boolean = true,
        sampleRate: Int = 24000,
        bitRate: Int = 24000,
    ) {
        encode(readAll(pcm), silk, sampleRate, bitRate, tencent)
    }

    private fun encode(
        pcm: ByteArray,
        silk: OutputStream,
        sampleRate: Int,
        bitRate: Int,
        tencent: Boolean,
    ) {
        validate(sampleRate, bitRate)
        writeHeader(silk, tencent)

        val encoder = SilkEncoder()
        encoder.configure(
            EncoderOptions(
                apiSampleRate = sampleRate,
                maxInternalSampleRate = 24000,
                packetSize = 20,
                bitRate = bitRate,
                packetLossPercentage = 0,
                complexity = 2,
                useInBandFec = false,
                useDtx = false,
            )
        )

        val apiPCM = decodePCM16LE(pcm)
        val internalRate = encoder.state.frameGeometry.samplingRateKHz * 1000
        val internalSampleCount = apiPCM.size.toLong() * internalRate / sampleRate
        val resampledPCM = resample(apiPCM, sampleRate, internalRate)
        check(internalSampleCount <= Int.MAX_VALUE && resampledPCM.size >= internalSampleCount)
        val internalPCM = if (resampledPCM.size == internalSampleCount.toInt()) {
            resampledPCM
        } else {
            resampledPCM.copyOf(internalSampleCount.toInt())
        }
        val frameLength = encoder.state.frameGeometry.frameLength
        val frame = IntArray(frameLength)
        val encoded = IntArray(MAX_ARITHM_BYTES)
        val encodedLength = RefInt(encoded.size)
        val payload = ByteArray(encoded.size)
        val payloadLengthBytes = ByteArray(2)
        var offset = 0
        while (offset + frameLength <= internalPCM.size) {
            internalPCM.copyInto(frame, 0, offset, offset + frameLength)
            encodeFrame(
                encoder.state,
                encoded,
                encodedLength,
                frame,
                payload,
                payloadLengthBytes,
                silk,
            )
            offset += frameLength
        }

        if (!tencent) silk.write(byteArrayOf(0xff.toByte(), 0xff.toByte()))
        silk.flush()
    }

    private fun writeHeader(silk: OutputStream, tencent: Boolean) {
        if (tencent) silk.write(0x02)
        silk.write(silkHeader)
    }

    private fun encodeFrame(
        encoder: EncoderState,
        encoded: IntArray,
        encodedLength: RefInt,
        frame: IntArray,
        payload: ByteArray,
        payloadLengthBytes: ByteArray,
        silk: OutputStream,
    ) {
        encodedLength.value = encoded.size
        val result = encodeFrame(encoder, encoded, encodedLength, frame)
        check(result == 0 && encodedLength.value in 1..encoded.size) {
            "SILK encoder returned $result with ${encodedLength.value} bytes"
        }

        val payloadLength = encodedLength.value
        for (index in 0 until payloadLength) payload[index] = encoded[index].toByte()
        payloadLengthBytes[0] = payloadLength.toByte()
        payloadLengthBytes[1] = (payloadLength ushr 8).toByte()
        silk.write(payloadLengthBytes)
        silk.write(payload, 0, payloadLength)
    }

    private fun readAll(source: InputStream): ByteArray {
        val output = ByteArrayOutputStream(64 * 1024)
        source.copyTo(output)
        return output.toByteArray()
    }

    private fun decodePCM16LE(bytes: ByteArray): IntArray {
        val samples = IntArray(bytes.size / 2)
        for (index in samples.indices) {
            val low = bytes[index * 2].toInt() and 0xff
            val high = bytes[index * 2 + 1].toInt()
            samples[index] = (high shl 8) or low
        }
        return samples
    }

    private fun validate(sampleRate: Int, bitRate: Int) {
        require(sampleRate in supportedSampleRates) {
            "Unsupported sample rate $sampleRate; expected one of ${supportedSampleRates.joinToString()}"
        }
        require(bitRate > 0) {
            "bitRate must be positive"
        }
    }
}
