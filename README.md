# silkt

silkt is a pure Kotlin SILK v3 encoder. Its encoder and
resampler core is a machine-assisted Kotlin port of
[tsilk](https://github.com/SnowLuma/tsilk), which is itself a TypeScript port
of the SILK fixed-point reference implementation.

## Usage

Input PCM must be mono, signed 16-bit, little endian. Supported input rates are
8000, 12000, 16000, and 24000 Hz. Packets are 20 ms. The output is a SILK v3 stream.

```kotlin
val silk = SilkCoder.encode(
    pcm = pcmBytes,
    sampleRate = 24000,
    bitRate = 24000,
    tencent = true,
)

SilkCoder.encode(pcmInputStream, silkOutputStream, sampleRate, bitRate)
```

## Development

```bash
./gradlew test
```
