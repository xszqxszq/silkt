package xyz.xszq.silkt.internal.codec

import xyz.xszq.silkt.internal.quantization.NoiseShapeQuantizerContext
import xyz.xszq.silkt.internal.signal.NoiseShapingContext
import xyz.xszq.silkt.internal.signal.PitchAnalysisContext
import xyz.xszq.silkt.internal.signal.PredictionAnalysisContext
import xyz.xszq.silkt.internal.signal.VoiceActivityState

internal class EncoderAnalysisState {
    val voiceActivity: VoiceActivityState = VoiceActivityState()
    val pitch: PitchAnalysisContext = PitchAnalysisContext()
    val noiseShaping: NoiseShapingContext = NoiseShapingContext()
    val prediction: PredictionAnalysisContext = PredictionAnalysisContext()
    val noiseShapeQuantization: NoiseShapeQuantizerContext = NoiseShapeQuantizerContext()
}
