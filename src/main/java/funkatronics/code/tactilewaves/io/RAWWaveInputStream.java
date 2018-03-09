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
 *  Class for acquiring audio from standard RAW Audio (*.raw) audio files.
 *  <p>
 *      Sets up a Java {@code InputStream} with the audio data and a {@link WaveFormat} with the
 *      audio format.
 *  </p>
 *
 * @author Marco Martinez
 */

public class RAWWaveInputStream extends WaveInputStream {

    private int mBlockAlign;

    private WaveFormat mFormat;
    private InputStream mStream;

    /**
     * Creates an {@code RAWWaveInputStream} Object to read audio data from a RAW Audio file
     * (*.raw)
     *
     * @param rawFile the File to read from
     * @param bigEndian true for Big Endian byte order, false for Little Endian
     * @param sampleRate the sample rate of the audio in the file
     * @param bitDepth the bit depth (bits/sample) of the audio in the file
     * @param numChannels the number of audio channels in the file
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public RAWWaveInputStream(File rawFile, boolean bigEndian, int sampleRate, int bitDepth, int numChannels) throws IOException{
        mFormat = new WaveFormat(WaveFormat.ENCODING_PCM_SIGNED, bigEndian, sampleRate, bitDepth, numChannels);
        mBlockAlign = numChannels*bitDepth/8;
        open(rawFile);
    }

    // Open the ile stream for reading
    private void open(File rawFile) throws IOException{
        this.mStream = new FileInputStream(rawFile);
    }

    public byte read() throws IOException{
        return (byte) mStream.read();
    }

    public float readSample() throws IOException{
        byte[] sample = new byte[mBlockAlign];
        int read = mStream.read(sample, 0, mBlockAlign);
        if(read == -1) return -2; // float samples are limited to [-1, 1] so -2 is used as EOF indicator
        float[] f = WaveFloatConverter.toMonoFloatArray(sample, mFormat);
        return f[0];
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
                "Raw Audio File %s", mFormat.toString());
    }
}
