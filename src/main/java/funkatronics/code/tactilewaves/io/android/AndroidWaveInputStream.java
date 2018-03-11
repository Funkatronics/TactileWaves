/*
 * ----------------------------------------------------------------------------)
 *                                                                            /
 *     _________  ________  ________ _________  ___  ___       _______       (
 *   |\___   ___\\   __  \|\   ____\\___   ___\\  \|\  \     |\  ___ \        \
 *   \|___ \  \_\ \  \|\  \ \  \___\|___ \  \_\ \  \ \  \    \ \   __/|        )
 *        \ \  \ \ \   __  \ \  \       \ \  \ \ \  \ \  \    \ \  \_|/__     /
 *         \ \  \ \ \  \ \  \ \  \____   \ \  \ \ \  \ \  \____\ \  \_|\ \   (
 *          \ \__\ \ \__\ \__\ \_______\  \ \__\ \ \__\ \_______\ \_______\   \
 *           \|__|  \|__|\|__|\|_______|   \|__|  \|__|\|_______|\|_______|    )
 *        ___       __   ________  ___      ___ _______   ________            /
 *       |\  \     |\  \|\   __  \|\  \    /  /|\  ___ \ |\   ____\          (
 *       \ \  \    \ \  \ \  \|\  \ \  \  /  / | \   __/|\ \  \___|_          \
 *        \ \  \  __\ \  \ \   __  \ \  \/  / / \ \  \_|/_\ \_____  \          )
 *         \ \  \|\__\_\  \ \  \ \  \ \    / /   \ \  \_|\ \|____|\  \        /
 *          \ \____________\ \__\ \__\ \__/ /     \ \_______\____\_\  \      (
 *           \|____________|\|__|\|__|\|__|/       \|_______|\_________\      \
 *                                                          \|_________|       )
 *                                                                            /
 *                                                                           (
 *                                                                            \
 * ----------------------------------------------------------------------------)
 *                                                                            /
 *                                                                           (
 *           Tactile Waves - a Java library for sensory substitution          \
 *          developed by Marco Martinez at Colorado State University.          )
 *                                                                            /
 *                                                                           (
 *                                                                            \
 * ----------------------------------------------------------------------------)
 *                                                                            /
 *                                                                           (
 *                    Copyright (C) 2017 Marco Martinez                       \
 *                                                                             )
 *                       funkatronicsmail@gmail.com                           /
 *                                                                           (
 *                                                                            \
 * ----------------------------------------------------------------------------)
 *                                                                            /
 *    This program is free software: you can redistribute it and/or modify   (
 *    it under the terms of the GNU General Public License as published by    \
 *     the Free Software Foundation, either version 3 of the License, or       )
 *                   (at your option) any later version.                      /
 *                                                                           (
 *       This program is distributed in the hope that it will be useful,      \
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of        )
 *        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the        /
 *                 GNU General Public License for more details               (
 *                                                                            \
 *     You should have received a copy of the GNU General Public License       )
 *    along with this program. If not, see <http://www.gnu.org/licenses/>.    /
 *                                                                           (
 *                                                                            \
 * ____________________________________________________________________________)
 */

package funkatronics.code.tactilewaves.io.android;

import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.os.Build;
import android.util.Log;

import java.io.IOException;

import funkatronics.code.tactilewaves.io.WaveFormat;
import funkatronics.code.tactilewaves.io.WaveInputStream;

/**
 *  Class for acquiring audio from the Android device's microphone.
 *  <p>
 *      Extends the {@link WaveInputStream} Abstract Class to wrap the Android AudioRecorder Object into
 *      a {@link WaveInputStream}.
 *  </p>
 *  <p>
 *      The highest sample rate available will be used, up to 44.1 kHz (The default Android sample
 *      rate).
 *  </p>
 *  <p>
 *      If the Android device being used has Android API level 21 (Lollipop) or higher, the class
 *      will attempt to use ENCODING_PCM_FLOAT to read floats directly from the OS. Otherwise, the
 *      Android default ENCODING_PCM_16BIT will be used and float conversion will happen after
 *      reading.
 *  </p>
 *
 * @author Marco Martinez
 *
 * @see WaveInputStream
 */
public class AndroidWaveInputStream extends WaveInputStream {

    private static String TAG = AndroidWaveInputStream.class.getSimpleName();

    // The stream to record audio from
    private AudioRecord mStream;
    // The format o the audio being acquired
    private WaveFormat mFormat;

    // The size of the audio buffer being read from, in samples
    private int mBufferSizeInSamples;

    /**
     * Creates an {@code AndroidWaveInputStream} Object and attempts to setup an {@code AudioRecord} object
     * to pipe audio data from the devices microphone
     *
     * @param bufferSizeInSamples the length of the read buffer, in samples
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public AndroidWaveInputStream(int bufferSizeInSamples) throws IOException{
        boolean result = getValidStream(bufferSizeInSamples);
        if(result) Log.i(TAG, "AndroidWaveStream initialized with " + mFormat.toString());
        else throw new IOException("AndroidWaveStream initialization FAILED");
        //start();
    }

    // Gets a valid input stream to record audio from the device hardware
    private boolean getValidStream(int bufferSizeInSamples) throws IOException{
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.LOLLIPOP) {
            for (int rate : new int[] {44100, 22050, 16000, 11025, 8000}) {
                int minBufferSizeInBytes = AudioRecord.getMinBufferSize(rate, AudioFormat.CHANNEL_IN_MONO, AudioFormat.ENCODING_PCM_FLOAT);
                if (minBufferSizeInBytes > 0) {
                    if(minBufferSizeInBytes/4 > bufferSizeInSamples) throw new IOException("buffer size must be at least " + minBufferSizeInBytes/4 + " samples long");
                    // buffer size is valid, Sample rate supported
                    Log.i(TAG, "sampleRate =  " + rate);
                    Log.i(TAG, "bufferSizeInSamples =  " + bufferSizeInSamples);
                    mFormat = new WaveFormat(WaveFormat.ENCODING_PCM_SIGNED, true, rate, 32, 1);
                    this.mBufferSizeInSamples = bufferSizeInSamples;
                    mStream = new AudioRecord(MediaRecorder.AudioSource.MIC,
                            mFormat.getSampleRate(),
                            AudioFormat.CHANNEL_IN_MONO,
                            AudioFormat.ENCODING_PCM_FLOAT,
                            4*bufferSizeInSamples);
                    break;
                }
            }
        } else {
            for (int rate : new int[] {44100, 22050, 16000, 11025, 8000}) {
                int minBufferSizeInBytes = AudioRecord.getMinBufferSize(rate, AudioFormat.CHANNEL_IN_MONO, AudioFormat.ENCODING_PCM_16BIT);
                if (minBufferSizeInBytes > 0) {
                    if(minBufferSizeInBytes/2 > bufferSizeInSamples) throw new IOException("buffer size must be at least " + minBufferSizeInBytes/2 + " samples long");
                    // buffer size is valid, Sample rate supported
                    Log.i(TAG, "sampleRate =  " + rate);
                    Log.i(TAG, "bufferSizeInSmaples =  " + bufferSizeInSamples);
                    mFormat = new WaveFormat(WaveFormat.ENCODING_PCM_SIGNED, true, rate, 16, 1);
                    this.mBufferSizeInSamples = bufferSizeInSamples;
                    mStream = new AudioRecord(MediaRecorder.AudioSource.MIC,
                            mFormat.getSampleRate(),
                            AudioFormat.CHANNEL_IN_MONO,
                            AudioFormat.ENCODING_PCM_16BIT,
                            2*bufferSizeInSamples);
                    break;
                }
            }
        }
        if(mStream.getState() == AudioRecord.STATE_INITIALIZED)
            return true;
        else
            return false;
    }

    public byte read() {
        int read;
        byte[] b = new byte[1];
        read = mStream.read(b, 0, 1);
        if(read > 0) return b[0];
        return -1;
    }

    public float readSample() {
        int read;
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.M
                && mStream.getAudioFormat() == AudioFormat.ENCODING_PCM_FLOAT) {
            float[] f = new float[1];
            read = mStream.read(f, 0, 1, AudioRecord.READ_BLOCKING);
            if(read > 0) return f[0];
        } else if (mStream.getAudioFormat() == AudioFormat.ENCODING_PCM_16BIT) {
            short[] s = new short[1];
            read = mStream.read(s, 0, 1);
            if(read > 0) return (float) (s[0] > 0 ? s[0] / Short.MAX_VALUE : s[0] / (Short.MAX_VALUE+1));
        }
        return -2;
    }

    public WaveFormat getFormat() {
        return mFormat;
    }

    /**
     * Get the minimum buffer size that was reported by the OS for the chosen sampleRate
     *
     * @return the minimum buffer size reported by the OS
     */
    public int getBufferSizeInSamples() {
        return mBufferSizeInSamples;
    }

    /**
     * Starts recording from the {@code AudioRecord} object. This will begin recording audio from the
     * device microphone
     */
    @Override
    public void start() {
        mStream.startRecording();
    }

    /**
     * Stops recording from the {@code AudioRecord} object. This will stop recording from the device
     * microphone (if it was already started).
     */
    @Override
    public void stop() {
        mStream.stop();
    }

    public void close() {
        mStream.release();
        mStream = null;
    }
}
