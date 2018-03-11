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
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Locale;

import funkatronics.code.tactilewaves.dsp.WaveFloatConverter;

/**
 *  Class for acquiring audio from standard WAVE (*.wav) audio files.
 *  <p>
 *      Automatically strips the WAVE file headers, reads the format parameters, and sets up a
 *      Java {@code InputStream} with the audio data and a {@link WaveFormat} with the audio
 *      format.
 *  </p>
 *
 */

public class MSWAVEInputStream extends WaveInputStream{

    private static String TAG = MSWAVEInputStream.class.getSimpleName();

    // The number of bytes used for 1 sample (including all channels)
    private int mBlockAlign;

    // The format of the Wave file
    private WaveFormat mFormat;
    // The file input stream to read from
    private InputStream mStream;

    // The number of bytes available to be read
    private long mBytesAvailable;

    /**
     * Creates an {@code MSWAVEInputStream} Object to read audio data from a standard WAVE file
     * (*.wav)
     *
     * @param waveFile the .wav file to read audio from
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public MSWAVEInputStream(File waveFile) throws IOException{
        open(waveFile);
    }

    /**
     * Returns the number of bytes available to be read from the stream
     *
     * @return the number of bytes available
     */
    public long bytesAvailable() {
        return mBytesAvailable;
    }

    /**
     * Returns the number of samples available to be read from the stream
     *
     * @return the number of samples available
     */
    public long samplesAvailable() {
        return mBytesAvailable / mBlockAlign;
    }

    private void open(File waveFile) throws IOException{
        this.mStream = new FileInputStream(waveFile);
        readHeader();
    }

    private void readHeader() throws IOException {
        readRIFF();
        readFMT();
        readDATA();
    }

    // Write the RIFF/WAVE portion of the WAVE header
    private void readRIFF() throws IOException{
        String riff = readID();
        if(!riff.equals("RIFF")) throw new IOException("Wave file header invalid");

        mBytesAvailable = readInt();

        String wave = readID();
        if(!wave.equals("WAVE")) throw new IOException("File is not in WAVE format");
    }

    // Write the format portion of the WAVE header
    private void readFMT() throws IOException{
        String Subchunk1ID = readID();
        if(!Subchunk1ID.equals("fmt ")) throw new IOException("Wave file format invalid");

        int Subchunk1Size = readInt();

        int fmt = readShort();
        if(fmt != 1) throw new IOException("Wave file is compressed, cannot read");

        int numChannels = readShort();

        int sampleRate = readInt();

        int byteRate = readInt();

        mBlockAlign = readShort();

        int bitsPerSample = readShort();

        if(byteRate != (sampleRate * numChannels * bitsPerSample/8))
            throw new IOException("Wave file byte rate invalid: " + byteRate);
        if(mBlockAlign != (numChannels * bitsPerSample /8))
            throw new IOException("Wave file block align invalid: " + mBlockAlign);

        mFormat = WaveFormat.MSWAVEFormat(sampleRate, bitsPerSample, numChannels);
    }

    // Write the data portion of the WAVE header
    private void readDATA() throws IOException{
        String Subchunk2ID = readID();
        if(!Subchunk2ID.equals("data")) throw new IOException("Wave data invalid");

        int numBytes = readInt();


        if(numBytes != mBytesAvailable - 36)
            throw new IOException("Wave header reports inconsistent data size");

        mBytesAvailable = numBytes;
    }

    // Read an ID (4 char String) from the stream
    private String readID() throws IOException {
        byte[] bytes = new byte[4];
        mStream.read(bytes);
        return new String(bytes);
    }

    // Read a short (2 bytes) from the stream
    private short readShort() throws IOException {
        byte[] bytes = new byte[2];
        mStream.read(bytes);
        return (short)((bytes[0] & 0xff)
                    | ((bytes[1] & 0xff) << 8 ));
    }

    // Read an int (4 bytes) from the stream
    private int readInt() throws IOException {
        byte[] bytes = new byte[4];
        mStream.read(bytes);
        return     (bytes[0] & 0xff)
                | ((bytes[1] & 0xff) << 8 )
                | ((bytes[2] & 0xff) << 16 )
                | ((bytes[3] & 0xff) << 24 );
    }

    public byte read() throws IOException{
        byte b = (byte) mStream.read();
        mBytesAvailable--;
        return b;
    }

    public float readSample() throws IOException{
        byte[] sample = new byte[mBlockAlign];
        int read = mStream.read(sample, 0, mBlockAlign);
        if(read == -1) return -2; // float samples are limited to [-1, 1] so -2 is used as EOF indicator
        mBytesAvailable -= mBlockAlign;
        float f = WaveFloatConverter.toMonoFloat(sample, mFormat);
        return f;
    }

    public WaveFormat getFormat() {
        return mFormat;
    }

    public void close() throws IOException {
        mStream.close();
        mStream = null;
    }

    @Override
    public String toString() {
        return String.format(Locale.getDefault(),
                "Wave File with %d samples. %s",
                samplesAvailable(), mFormat.toString());
    }
}
