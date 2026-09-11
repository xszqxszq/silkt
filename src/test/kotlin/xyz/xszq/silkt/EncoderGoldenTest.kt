package xyz.xszq.silkt

import xyz.xszq.silkt.internal.codec.EncoderOptions
import xyz.xszq.silkt.internal.codec.SilkEncoder
import xyz.xszq.silkt.internal.codec.encodeFrame
import xyz.xszq.silkt.internal.signal.resample
import xyz.xszq.silkt.internal.util.RefInt
import kotlin.math.PI
import kotlin.math.sin
import kotlin.test.Test
import kotlin.test.assertEquals

class EncoderGoldenTest {
    @Test
    fun matchesTsilkSineWaveAt24k() {
        val apiPCM = IntArray(6 * 480) {
            (10000.0 * sin(2.0 * PI * 440.0 * it / 24000.0)).toInt()
        }
        val pcm = resample(apiPCM, 24000, 16000)
        assertEquals(1920, pcm.size)

        val encoder = SilkEncoder()
        encoder.configure(
            EncoderOptions(
                apiSampleRate = 24000,
                maxInternalSampleRate = 24000,
                packetSize = 20,
                bitRate = 24000,
                packetLossPercentage = 0,
                complexity = 2,
                useInBandFec = false,
                useDtx = false,
            )
        )

        val expectedHex = listOf(
            "8a3c5bd0d348d9f4288f59db703a57e7b132fdcf772ea6ed8d9b1d963800a67dbec605fd7c3f",
            "81b521553c5ccb2fb1ff8fc4a1be0b8aad3fbd4dfe72c517cdb6d8c27312e8ecfce514d04792e0bf",
            "8151595f39a1dd1383811bf66e6ee883e452a82121cdd6c1119475f5535c1800548b3600d6b045" +
                    "a201ceffe1f8731730110f0bd0b3c58ae00218f9b0982945e2ea30f976503e3f",
            "8151595f39a1dd13838007a047b9a9a810c9626ffb83895cf18d550992072c62725e5070bff93ede7d" +
                    "e1df7cd092d7ecb3b57189ff54a73357",
            "81553c17b3ea745ee7b578aa613e030e435af451f8e7b6675a67efb97137923a1b668df1dc0e5d67d83caa615e8062e5d127be9f",
            "815159a6b9d766f8bd430f67d8a1a5a24139e2dfeb4c837c6d8708b5f6457cfe7f68f72c73f6db24" +
                    "cb7f987c7068b0436bdb8013f35816abb5bfc0d045cb",
        )
        val payloads = mutableListOf<IntArray>()
        var offset = 0
        while (offset < pcm.size) {
            val output = IntArray(1250)
            val byteCount = RefInt(output.size)
            encodeFrame(encoder.state, output, byteCount, pcm.copyOfRange(offset, offset + 320))
            payloads += output.copyOf(byteCount.value)
            offset += 320
        }

        val expected = expectedHex.joinToString("")
        val actual = payloads.joinToString("") { frame ->
            frame.joinToString("") { it.toString(16).padStart(2, '0') }
        }
        assertEquals(expected, actual)
    }

}
