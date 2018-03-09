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

package funkatronics.code.tactilewaves.dsp;

import java.util.Arrays;

import funkatronics.code.tactilewaves.io.WaveFormat;

/**
 * Class to convert audio encoded in bytes to floats and vice versa.
 *
 * @author Marco Martinez
 *
 * @see funkatronics.code.tactilewaves.io.MSWAVEInputStream
 * @see funkatronics.code.tactilewaves.io.RAWWaveInputStream
 * @see funkatronics.code.tactilewaves.io.WaveWriter
 */

public class WaveFloatConverter {

    /**
     * Convert an array of bytes to a single sample.
     * <p>
     *     The length of the array of bytes should equal the bytes per sample * the number of
     *     channels of the {@link WaveFormat}. If it is less than this length, an exception is
     *     thrown. If it is longer, the extra bytes will be ignored.
     * </p>
     *
     * @param b the array of bytes representing a single sample
     * @param format the {@code WaveFormat} of the audio bytes
     *
     * @return a sample (single precision)
     */
    public static float toMonoFloat(byte[] b, WaveFormat format){
        if(format.isBigEndian())
            return toMonoFloatBigEndian(b, format.getChannels(), format.getBitDepth());
        else
            return toMonoFloatLittleEndian(b, format.getChannels(), format.getBitDepth());
    }

    /**
     * Convert an array of bytes to a single sample.
     * <p>
     *     The length of the array of bytes should equal the bytes per sample * the number of
     *     channels of the {@link WaveFormat}. If it is less than this length, an exception is
     *     thrown. If it is longer, the extra bytes will be ignored.
     * </p>
     *
     * @param b the array of bytes representing a single sample
     * @param format the {@code WaveFormat} of the audio bytes
     *
     * @return a sample (double precision)
     */
    public static double toMonoDouble(byte[] b, WaveFormat format){
        if(format.isBigEndian())
            return toMonoDoubleBigEndian(b, format.getChannels(), format.getBitDepth());
        else
            return toMonoDoubleLittleEndian(b, format.getChannels(), format.getBitDepth());
    }

    /**
     * Convert a single sample to an array of bytes.
     *
     * @param f the sample to convert
     * @param format the {@code WaveFormat} of the audio sample
     *
     * @return an array of bytes representing the sample according to the {@code WaveFormat}
     */
    public static byte[] toBytes(float f, WaveFormat format) {
        if(format.isBigEndian())
            return toBytesBigEndian(f, format.getChannels(), format.getBitDepth());
        else
            return toBytesLittleEndian(f, format.getChannels(), format.getBitDepth());
    }

    /**
     * Convert a single sample to an array of bytes.
     *
     * @param d the sample to convert
     * @param format the {@code WaveFormat} of the audio sample
     *
     * @return an array of bytes representing the sample according to the {@code WaveFormat}
     */
    public static byte[] toBytes(double d, WaveFormat format) {
        if(format.isBigEndian())
            return toBytesBigEndian(d, format.getChannels(), format.getBitDepth());
        else
            return toBytesLittleEndian(d, format.getChannels(), format.getBitDepth());
    }

    /**
     * Convert an array of bytes to an array of samples.
     * <p>
     *     The length of the array of bytes should be divisible by the bytes per sample * the number
     *     of channels of the {@link WaveFormat}. If it is not the extra bytes will be ignored.
     * </p>
     *
     * @param buffer the array of bytes representing a bunch of samples
     * @param format the {@code WaveFormat} of the audio bytes
     *
     * @return an array of samples (single precision)
     */
    public static float[] toMonoFloatArray(byte[] buffer, WaveFormat format) {
        if(format.isBigEndian())
            return toMonoFloatArrayBigEndian(buffer, format.getChannels(), format.getBitDepth());
        else
            return toMonoFloatArrayLittleEndian(buffer, format.getChannels(), format.getBitDepth());
    }

    /**
     * Convert an array of bytes to an array of samples.
     * <p>
     *     The length of the array of bytes should be divisible by the bytes per sample * the number
     *     of channels of the {@link WaveFormat}. If it is not the extra bytes will be ignored.
     * </p>
     *
     * @param buffer the array of bytes representing a bunch of samples
     * @param format the {@code WaveFormat} of the audio bytes
     *
     * @return an array of samples (double precision)
     */
    public static double[] toMonoDoubleArray(byte[] buffer, WaveFormat format) {
        if(format.isBigEndian())
            return toMonoDoubleArrayBigEndian(buffer, format.getChannels(), format.getBitDepth());
        else
            return toMonoDoubleArrayLittleEndian(buffer, format.getChannels(), format.getBitDepth());

    }

    // convert an array of bytes to a float
    private static float toMonoFloatLittleEndian(byte[] b, int channels, int bitDepth){
        int blockAlign = channels*bitDepth/8;
        if(b.length < blockAlign)
            throw new IllegalArgumentException("The length of the byte buffer should be equal to channels*bitDepth/8");
        float sample = 0;
        long temp = 0;
        for(int c = 0; c < channels; c++) {
            int shift = c * (blockAlign / channels);
            for (int i = 0; i < blockAlign / channels; i++)
                temp = temp | ((b[i + shift] & 0xffL) << 8L * i);
            sample += toFloat(temp, bitDepth);
        }
        sample /= channels;
        return sample;
    }

    // convert an array of bytes to a float
    private static double toMonoDoubleLittleEndian(byte[] b, int channels, int bitDepth){
        int blockAlign = channels*bitDepth/8;
        if(b.length < blockAlign)
            throw new IllegalArgumentException("The length of the byte buffer should be equal to channels*bitDepth/8");
        double sample = 0;
        long temp = 0;
        for(int c = 0; c < channels; c++) {
            int shift = c * (blockAlign/channels);
            for (int i = 0; i < blockAlign/channels; i++)
                temp = temp | ((b[i + shift] & 0xffL) << 8L * i);
            sample += toDouble(temp, bitDepth);
        }
        sample /= channels;
        return sample;
    }

    // convert a float to an array of bytes
    private static byte[] toBytesLittleEndian(float f, int channels, int bitDepth) {
        long sampleAsLong = toLong(f, bitDepth);
        int blockAlign = channels*(bitDepth/8);
        byte[] bytes = new byte[blockAlign];
        for(int c = 0; c < channels; c++) {
            int shift = c*(blockAlign/channels);
            for (int b = 0; b < blockAlign/channels; b++) {
                bytes[b + shift] = (byte) ((sampleAsLong >> 8L * b) & 0xff);
            }
        }
        return bytes;
    }

    // convert a float to an array of bytes
    private static byte[] toBytesLittleEndian(double d, int channels, int bitDepth) {
        long sampleAsLong = toLong(d, bitDepth);
        int blockAlign = channels*(bitDepth/8);
        byte[] bytes = new byte[blockAlign];
        for(int c = 0; c < channels; c++) {
            int shift = c*(blockAlign/channels);
            for (int b = 0; b < blockAlign/channels; b++)
                bytes[b + shift] = (byte) ((sampleAsLong >> 8L * b) & 0xffL);
        }
        return bytes;
    }

    // convert an array of bytes to a float
    private static float toMonoFloatBigEndian(byte[] b, int channels, int bitDepth){
        int blockAlign = channels*bitDepth/8;
        if(b.length < blockAlign)
            throw new IllegalArgumentException("The length of the byte buffer should be equal to channels*bitDepth/8");
        float sample = 0;
        long temp = 0;
        for(int c = 0; c < channels; c++) {
            int shift = c * (blockAlign / channels);
            for (int i = 0; i < blockAlign / channels; i++)
                temp = temp | ((b[blockAlign - (i + shift) - 1] & 0xffL) << 8L * i);
            sample += toFloat(temp, bitDepth);
        }
        sample /= channels;
        return sample;
    }

    // convert an array of bytes to a float
    private static double toMonoDoubleBigEndian(byte[] b, int channels, int bitDepth){
        int blockAlign = channels*bitDepth/8;
        if(b.length < blockAlign)
            throw new IllegalArgumentException("The length of the byte buffer should be equal to channels*bitDepth/8");
        double sample = 0;
        long temp = 0;
        for(int c = 0; c < channels; c++) {
            int shift = c * (blockAlign / channels);
            for (int i = 0; i < blockAlign / channels; i++)
                temp = temp | ((b[blockAlign - (i + shift) - 1] & 0xffL) << 8L * i);
            sample += toDouble(temp, bitDepth);
        }
        sample /= channels;
        return sample;
    }

    // convert a float to an array of bytes
    private static byte[] toBytesBigEndian(float f, int channels, int bitDepth) {
        long sampleAsLong = toLong(f, bitDepth);
        int blockAlign = channels*(bitDepth/8);
        byte[] bytes = new byte[blockAlign];
        for(int c = 0; c < channels; c++) {
            int shift = c*(blockAlign/channels);
            for (int b = 0; b < blockAlign/channels; b++) {
                bytes[blockAlign - (b + shift) - 1] = (byte) ((sampleAsLong >> 8L * b) & 0xff);
            }
        }
        return bytes;
    }

    // convert a float to an array of bytes
    private static byte[] toBytesBigEndian(double d, int channels, int bitDepth) {
        long sampleAsLong = toLong(d, bitDepth);
        int blockAlign = channels*(bitDepth/8);
        byte[] bytes = new byte[blockAlign];
        for(int c = 0; c < channels; c++) {
            int shift = c*(blockAlign/channels);
            for (int b = 0; b < blockAlign/channels; b++) {
                bytes[blockAlign - (b + shift) - 1] = (byte) ((sampleAsLong >> 8L * b) & 0xff);
            }
        }
        return bytes;
    }

    // convert a byte array to an array of floats
    private static float[] toMonoFloatArrayLittleEndian(byte[] buffer, int channels, int bitDepth) {
        int blockAlign = channels*bitDepth/8;
        float[] samples = new float[buffer.length/blockAlign];
        for(int i = 0; i < samples.length; i++) {
            samples[i] = toMonoFloatLittleEndian(
                    Arrays.copyOfRange(buffer,i*blockAlign, blockAlign), channels, bitDepth);
        }
        return samples;
    }

    // convert a byte array to an array of doubles
    private static double[] toMonoDoubleArrayLittleEndian(byte[] buffer, int channels, int bitDepth) {
        int blockAlign = channels*bitDepth/8;
        double[] samples = new double[buffer.length/blockAlign];
        for(int i = 0; i < samples.length; i++) {
            samples[i] = toMonoDoubleLittleEndian(
                    Arrays.copyOfRange(buffer,i*blockAlign, blockAlign), channels, bitDepth);
        }
        return samples;
    }

    // convert a byte array to an array of floats
    private static float[] toMonoFloatArrayBigEndian(byte[] buffer, int channels, int bitDepth) {
        int blockAlign = channels*bitDepth/8;
        float[] samples = new float[buffer.length/blockAlign];
        for(int i = 0; i < samples.length; i++) {
            samples[i] = toMonoFloatBigEndian(
                    Arrays.copyOfRange(buffer,i*blockAlign, blockAlign), channels, bitDepth);
        }
        return samples;
    }

    // convert a byte array to an array of floats
    private static double[] toMonoDoubleArrayBigEndian(byte[] buffer, int channels, int bitDepth) {
        int blockAlign = channels*bitDepth/8;
        double[] samples = new double[buffer.length/blockAlign];
        for(int i = 0; i < samples.length; i++) {
            samples[i] = toMonoDoubleBigEndian(
                    Arrays.copyOfRange(buffer,i*blockAlign, blockAlign), channels, bitDepth);
        }
        return samples;
    }

    // convert a Long to a float, extending the sign
    private static float toFloat(long L, int bitDepth) {
        int extensionBits = Long.SIZE - bitDepth;
        L = (L << extensionBits) >> extensionBits;
        return (float) (L > 0 ? L / (Math.pow(2, bitDepth - 1)-1)
                : L / (Math.pow(2, bitDepth - 1)));
    }

    // convert a Long to a double, extending the sign
    private static double toDouble(long L, int bitDepth) {
        int extensionBits = Long.SIZE - bitDepth;
        L = (L << extensionBits) >> extensionBits;
        return L > 0 ? L / (Math.pow(2, bitDepth - 1)-1)
                : L / (Math.pow(2, bitDepth - 1));
    }

    // convert a float to a long,
    private static Long toLong(float f, int bitDepth) {
        long L = (long) (f * Math.pow(2, bitDepth - 1));
        return L;
    }

    //convert a double to a long
    private static Long toLong(double d, int bitDepth) {
        long L = (long) (d * Math.pow(2, bitDepth - 1));
        return L;
    }

}
