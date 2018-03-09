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

import java.io.IOException;

import funkatronics.code.tactilewaves.io.android.AndroidWaveInputStream;

/**
 * This abstract class is the superclass of all classes representing an audio input stream.
 * <p>
 *     Any object that is a subclass of InputStream must always provide methods that
 *     return the next byte of input, the next sample of input, and the {@link WaveFormat} in use.
 *     Additionally, subclasses must implement a method to close the underlying stream.
 * </p>
 *
 * @author Marco Martinez
 *
 * @see AndroidWaveInputStream
 * @see MSWAVEInputStream
 */

public abstract class WaveInputStream {

    /**
     * Reads the next byte from the input stream and returns it, or returns -1 if a byte could not
     * be read
     *
     * @return the next sample value in the stream or -2
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public abstract byte read() throws IOException;

    /**
     * Reads the next sample(float) from the input stream and returns it, or returns -2 if a sample
     * could not be read. If the stream contains multi-channel audio, the next sample from each
     * channel is read, and the average is returned.
     *
     * @return the next sample value in the stream or -2
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public abstract float readSample() throws IOException;

    /**
     * Reads floats of audio data from the input stream up to the length of the buffer
     *
     * @param f the byte buffer array to which data is read
     *
     * @return the number of bytes successfully read
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public int read(float[] f) throws IOException {
        return read(f, 0, f.length);
    }

    /**
     * Reads floats of audio data from the input stream up to the length of the buffer
     *
     * @param f the byte buffer array to which data is read
     * @param off location to begin writing to the array
     * @param len how many bytes to read (len > (b.length-off) ? len = b.length-off)
     *
     * @return the number of bytes successfully read
     *
     * @throws IOException if an input/output error occurs on the stream
     */
    public int read(float[] f, int off, int len) throws IOException {
        if(len > f.length-off) len = f.length - off;
        int read = 0;
        while(read < len) {
            float temp = readSample();
            if(temp < -1) return -1;
            f[off + read++] = temp;
        }
        return read;
    }

    /**
     * Get the {@link WaveFormat} associated with the underlying audio stream
     *
     * @return the {@link WaveFormat} associated with this audio stream
     */
    public abstract WaveFormat getFormat();

    /**
     * Start the underlying stream
     */
    public void start() {}

    /**
     * Stop the underlying stream
     */
    public void stop() {}

    /**
     * Close the underlying stream and release its resources
     *
     * @throws IOException if an input/output error occurs when closing the stream
     */
    public abstract void close() throws IOException;
}
