package xyz.xszq.silkt

import java.io.ByteArrayInputStream
import java.io.ByteArrayOutputStream
import kotlin.math.PI
import kotlin.math.sin
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class SilkCoderTest {
    @Test
    fun writesTencentContainerMatchingPayloadGolden() {
        val pcm = sinePCM(6)
        val silk = SilkCoder.encode(pcm)

        assertEquals(0x02, silk[0].toInt() and 0xff)
        assertEquals("#!SILK_V3", silk.copyOfRange(1, 10).decodeToString())
        assertEquals(
            listOf(38, 40, 71, 57, 52, 62),
            parseFrameLengths(silk),
        )
        assertTrue(silk.copyOfRange(silk.size - 2, silk.size).contentEquals(byteArrayOf(-1, -1)).not())
    }

    @Test
    fun streamAndByteArrayEncodersAgree() {
        val pcm = sinePCM(4)
        val streamed = ByteArrayOutputStream()
        SilkCoder.encode(ByteArrayInputStream(pcm), streamed)

        val expected = SilkCoder.encode(pcm)
        val actual = streamed.toByteArray()
        assertTrue(expected.contentEquals(actual))
    }

    @Test
    fun streamEncoderWritesStandardContainer() {
        val standard = ByteArrayOutputStream()
        SilkCoder.encode(ByteArrayInputStream(sinePCM(1)), standard, tencent = false)
        val silk = standard.toByteArray()

        assertEquals("#!SILK_V3", silk.copyOfRange(0, 9).decodeToString())
        assertTrue(silk.copyOfRange(silk.size - 2, silk.size).contentEquals(byteArrayOf(-1, -1)))
    }

    @Test
    fun writesStandardEndMarkerAndRejectsUnsupportedRates() {
        val silk = SilkCoder.encode(sinePCM(1), tencent = false)
        assertEquals("#!SILK_V3", silk.copyOfRange(0, 9).decodeToString())
        assertTrue(silk.copyOfRange(silk.size - 2, silk.size).contentEquals(byteArrayOf(-1, -1)))
        assertFailsWith<IllegalArgumentException> {
            SilkCoder.encode(ByteArray(0), sampleRate = 44100)
        }
    }

    @Test
    fun ignoresTrailingPartialFrameLikeReferenceEncoder() {
        val fullFrames = sinePCM(3)
        val partial = fullFrames + ByteArray(17)

        assertEquals(3, parseFrameLengths(SilkCoder.encode(partial)).size)
    }

    @Test
    fun resamplerFilterTailDoesNotCreateAnExtraFrame() {
        val pcm = sinePCM(40).copyOf(40 * 480 * 2 - 2)

        assertEquals(39, parseFrameLengths(SilkCoder.encode(pcm)).size)
    }

    @Test
    fun encodesAt24kHzInternalSampleRate() {
        val lengths = parseFrameLengths(SilkCoder.encode(sinePCM(3), bitRate = 50_000))

        assertEquals(3, lengths.size)
    }

    private fun sinePCM(frames: Int): ByteArray {
        val sampleCount = frames * 480
        return ByteArray(sampleCount * 2) { index ->
            val sample = (10000.0 * sin(2.0 * PI * 440.0 * (index / 2) / 24000.0)).toInt()
            if (index % 2 == 0) sample.toByte() else (sample shr 8).toByte()
        }
    }

    private fun parseFrameLengths(silk: ByteArray): List<Int> {
        val lengths = mutableListOf<Int>()
        var offset = 10
        while (offset + 2 <= silk.size) {
            val length = (silk[offset].toInt() and 0xff) or
                    ((silk[offset + 1].toInt() and 0xff) shl 8)
            if (length == 0xffff) break
            lengths += length
            offset += 2 + length
        }
        return lengths
    }
}
