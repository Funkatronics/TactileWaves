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

package funkatronics.code.tactilewaves.io;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Locale;

import funkatronics.code.tactilewaves.dsp.WaveFloatConverter;

/**
 * Class to write audio data to a Wave (.wav) file according to the format described
 * <a href="http://soundfile.sapp.org/doc/WaveFormat/">here</a>.
 *
 * @author Marco Martinez
 *
 * @see WaveFormat
 * @see MSWAVEInputStream
 */

public class WaveWriter {

    // The length (in bytes) of the file header
    private  static final int HEADER_LENGTH = 44;

    // Constant that represents the PCM encoding format of Wave files
    private static final short AUDIOFORMAT_PCM = 1;
    // Constant that represents the ALAW encoding format of Wave files
    private static final short AUDIOFORMAT_ALAW = 6;
    // Constant that represents the ULAW encoding format of Wave files
    private static final short AUDIOFORMAT_ULAW = 7;

    // The number of distinct channels in the written file
    private int mNumChannels;
    // The sample rate of the audio data written to the file
    private int mSampleRate;
    // The number of bits to use for each sample of audio data written to the file
    private int mBitDepth;

    // The number of samples written to the file - needed in the file header
    private int mNumSamples;

    // The format of the audio being written to the file
    private WaveFormat mFormat;

    // The file output
    private RandomAccessFile mOut;

    /**
     * Construct a new WaveWriter from a file name and {@link WaveFormat}
     *
     * @param waveFileName the name of the resulting .wav file (including the path)
     * @param format the format of audio to be written
     *
     * @throws IOException if an error occurs when opening the file
     */
    public WaveWriter(String waveFileName, WaveFormat format) throws IOException {
        this(new File(waveFileName), format.getChannels(), format.getSampleRate(), format.getBitDepth());
    }

    /**
     * Construct a new WaveWriter from a File and {@link WaveFormat}
     *
     * @param waveFile the File to write data to
     * @param format the format of audio to be written
     *
     * @throws IOException if an error occurs when opening the file
     */
    public WaveWriter(File waveFile, WaveFormat format) throws IOException {
        this(waveFile, format.getChannels(), format.getSampleRate(), format.getBitDepth());
    }

    /**
     * Construct a new WaveWriter from a File and the audio format parameters
     *
     * @param waveFile the File to write data to
     * @param numChannels the number of channels in the audio
     * @param sampleRate the sample rate (Hz) of the audio
     * @param bitsPerSample how many bits to use to store each sample
     *
     * @throws IOException if an error occurs when opening the file
     */
    public WaveWriter(File waveFile, int numChannels, int sampleRate, int bitsPerSample) throws IOException {
        mNumChannels = numChannels;
        mSampleRate = sampleRate;
        mBitDepth = bitsPerSample;
        mNumSamples = 0;
        mFormat = new WaveFormat(WaveFormat.ENCODING_PCM_SIGNED, false, sampleRate,
                bitsPerSample, numChannels);
        open(waveFile);
    }

    // Open the file output for writing
    private void open(File file) throws IOException{
        mOut = new RandomAccessFile(file, "rwd");
        mOut.seek(0);
        mOut.write(new byte[HEADER_LENGTH]);
    }

    /**
     * Writes the proper WAVE header to the file and closes it.
     *
     * @throws IOException if an error occurs when closing the file
     */
    public void close() throws IOException {
        mOut.seek(0);
        writeHeader();
        mOut.close();
        mOut = null;
    }

    // Write the WAVE header
    private void writeHeader() throws IOException {
        writeRIFF();
        writeFMT();
        writeDATA();
    }

    // Write the RIFF/WAVE portion of the WAVE header
    private void writeRIFF() throws IOException{
        int chunkSize = 36 + mNumSamples * mNumChannels * (mBitDepth/8);
        write("RIFF");
        write(chunkSize);
        write("WAVE");
    }

    // Write the format portion of the WAVE header
    private void writeFMT() throws IOException{
        int Subchunk1Size = 16; // 16 for PCM
        write("fmt ");
        write(Subchunk1Size);
        write((short)AUDIOFORMAT_PCM);
        write((short) mNumChannels);
        write(mSampleRate);
        write(mSampleRate * mNumChannels * mBitDepth/8);
        write((short) (mNumChannels * mBitDepth/8));
        write((short) mBitDepth);
    }

    // Write the data portion of the WAVE header
    private void writeDATA() throws IOException{
        int Subchunk2Size = mNumSamples * mNumChannels * (mBitDepth/8);
        write("data");
        write(Subchunk2Size);
    }

    /**
     * Write a single sample to the file
     *
     * @param sample the sample to write
     *
     * @throws IOException if an error occurs when writing to the file
     */
    public void writeSample(float sample) throws IOException {
        mOut.write(WaveFloatConverter.toBytes(sample, mFormat));
        mNumSamples++;
    }

    /**
     * Write a single sample to the file
     *
     * @param sample the sample to write
     *
     * @throws IOException if an error occurs when writing to the file
     */
    public void writeSample(double sample) throws IOException {
        mOut.write(WaveFloatConverter.toBytes(sample, mFormat));
        mNumSamples++;
    }

    /**
     * Write an array of samples to the file
     *
     * @param buffer the sample buffer/array to write
     *
     * @throws IOException if an error occurs when writing to the file
     */
    public void write(float[] buffer) throws IOException {
        for(float sample : buffer) {
            writeSample(sample);
        }
    }

    /**
     * Write an array of samples to the file
     *
     * @param buffer the sample buffer/array to write
     *
     * @throws IOException if an error occurs when writing to the file
     */
    public void write(double[] buffer) throws IOException {
        for(double sample : buffer) {
            writeSample(sample);
        }
    }

    // Write a string to the file - used to write the header
    private void write(String s) throws IOException {
        mOut.write(s.getBytes("UTF-8"), 0, s.length());
    }

    // Write an int to the file - used to write the header
    private void write(int i) throws IOException {
        mOut.write( i        & 0x00ff);
        mOut.write((i >>  8) & 0x00ff);
        mOut.write((i >> 16) & 0x00ff);
        mOut.write((i >> 24) & 0x00ff);
    }

    // Write a short to the file - used to write the header
    private void write(short s) throws IOException {
        mOut.write( s       & 0x00ff);
        mOut.write((s >> 8) & 0x00ff);
    }

    @Override
    public String toString() {
        return String.format(Locale.getDefault(),
                "Wave File with %d samples. %s",
                mNumSamples, mFormat.toString());
    }

}
