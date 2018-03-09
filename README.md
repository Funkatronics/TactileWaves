# Tactile Waves

Tactile Waves is an open source digital signal processing (DSP) library for sound-to-touch sensory 
substitution research and development. It has been designed to enable mobile devices to be used as 
the 'back end' of a sound-to-touch sensory substitution system. Using Tactile Waves, audio can be 
acquired and preprocessed, and important features can be extracted and sent to sensory substitution 
hardware worn on the body via Bluetooth. 

Tactile Waves can be useful even if you are'nt interested in sensory substitution. As required by 
sensory substitution applications, Tactile Waves provides a plethora of speech processing functions 
and audio stream management. These can be used to build audio analysis and feature extraction 
applications for Android and Java such as a Pitch Detector or Formant Frequency Estimator. Tactile
Waves could also be used as a preprocessor and feature extractor for an Automatic Speech Recognition
(ASR) system. The possibilities are endless, and it is my hope that users will be encouraged to push
the boundaries of the library (and contribute their extensions!).

## Installing Tactile Waves

Tactile Waves is available as a remote Android Library, as well as a standard JAR file.  

If you are developing an Android Application, the easiest way to use Tactile waves is to add a 
dependency to the remote library.

Add the following dependency to your apps build.gradle file:

```
dependencies {
    compile 'com.funkatronics.code:tactilewaves:1.0'
}
```

For a standard java application, simply download the provided JAR file and add it to your projects /libs folder.

## Getting Started with Tactile Waves

Tactile Waves provides a system for handling audio recording, buffering, windowing and processing in 
a few easy to use objects within the `tactilewaves.dsp` package. 

### WaveManager

The `WaveManager` object represents the core of Tactile Waves' audio engine. It reads audio samples 
from a `WaveInputStream`, and generates a `WaveFrame` for each frame of audio acquired. Each 
`WaveFrame` is passed through a chain of `WaveProcessor` objects that perform some useful processing
on the frame, and optionally store extracted audio features in the frame's feature set. Processed 
frames are then sent to any `WaveListener` objects that have been added to the `WaveManager`. These 
listeners can be used to trigger UI updates, or send extracted audio features to sensory 
substitution hardware.

### WaveInputStream

'WaveInputStream' is an abstract class that contains an underlying `InputStream` to acquire audio 
bytes from a stream and convert them audio samples. A `WaveManager` object is used to read a 
`WaveInputStream` and generate windowed frames of audio.

### WaveFrame

The `WaveFrame` object represents a single frame of audio. A `WaveManager` will create `WaveFrame` 
objects according the the buffer and overlap length provided. Each `WaveFrame` is then passed 
through the `WaveManager`'s processing chain

### WaveProcessor

Tactile Waves uses a `WaveProcessor` interface to build objects that can be inserted into a 
`WaveManager`'s processing chain. Any object that implements this interface can be added to the 
chain to perform some type of processing or analysis on each frame of audio.

### WaveFrameListener

After a `WaveFrame` has been passed through the entire processing chain, it is sent to any 
subscribed `WaveFrameListener` objects. This can be used to trigger events that are dependent on 
processed audio frames such as UI updates or Bluetooth transmission. 

### Putting It All Together
 
A `WaveManager` is initialized with a `WaveInputStream` to read audio from, along with a buffer and overlap length.

To obtain a stream from an Android device's microphone use:
```
WaveInputStream inputStream = new AndroidWaveInputStream(bufferSize);
```
A `WaveManager` can then be instantiated:
```
WaveManager waveManager = new WaveManager(inputStream, bufferSize, overlap);
```
A newly instantiated `WaveManager` will have an empty processing chain. To add processing to the 
chain, use the `addEffectToChain` method to add any object that implements the `WaveProcessor` 
interface. Tactile Waves provides several ready to use processors, and users can easily build their 
own.

For example, to estimate the pitch of incoming audio, a `PitchProcessor` can be used. Add a new 
`PitchProcessor` the the Wave Manager's processing chain:
```
waveManager.addEffectToChain(new PitchProcessor(bufferSize));
```
Add any listeners to the `WaveManager` using the `addListener()` method:
```
waveManager.addListener(listener);
```
The `WaveManager` is now ready to be started. Calling `start()` will start processing audio in a new
thread. Calling `stop()` will stop this thread.
```
waveManager.start();
```
`WaveManager` implements the `Runnable` interface allowing for more specialized thread management:
```
Thread thread = new Thread(waveManager);
thread.start();
// Do something with the thread
```

## Deployment

Tactile Waves has been tested on a Moto G5 Plus running Android 7.0 on a Snapdragon 625 as well as a
Windows laptop with an Intel i7-8550U processor. Real-time operation has been confirmed on these 
devices. More Android device tests will be performed in the near future, but any modern Android 
device should be more than capable of processing audio in real-time with Tactile Waves. 

The built-in `StopWatch` class can be used in a `WaveFrameListener` object to measure the elapsed
time between completed audio frames, to ensure real-time operation on your device/system.

## Built With

* [Android Studio](https://developer.android.com/studio/index.html)
* [JFrog Bintray](https://bintray.com/)
* [Maven](https://maven.apache.org/)
* [BlueCove](http://bluecove.org/)

## Contributing

Contributions to Tactile Waves are welcomed. If you would live to contribute to the code simply fork the project on GitHub and send me a pull request when you have code ready to be be included. 

## Versioning

Version 1.0 - 2018-03-09
* First public release


## Authors

* **Marco Martinez**(*Funkatronics*)

## License

This software is licensed under the GNU General Public License - see the[license.txt](license.txt)file for details

## Funding

This project was funded by a Creative Works Award from[CSU Ventures](http://csuventures.org/).

## Acknowledgments

* Joren Six for his[TarsosDSP](https://github.com/JorenSix/TarsosDSP)Library that provided an excellent learning resource in the early days of my thesis
* JJ Moritz at[Sapien, LLC](http://www.sapienllc.com/)for bringing me into CSU's sensory substitution research and providing assistance and resources throughout the project
* Dr. Leslie Stone-Roy and Dr. John Williams at Colorado State University for their advice and guidance
