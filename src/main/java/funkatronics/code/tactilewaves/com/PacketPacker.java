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

package funkatronics.code.tactilewaves.com;

import java.util.Arrays;
import java.util.BitSet;

import funkatronics.code.tactilewaves.com.android.BluetoothAndroid;

/**
 * Class to pack data into a data packet for Bluetooth transmission
 *
 * @author Marco Martinez
 *
 * @see BluetoothAndroid
 */

public class PacketPacker {

    private int mBitDepth;

    private double mMinValue;
    private double mRange;

    /**
     * Construct a new {@code PacketPacker}
     *
     * @param bitsPerChannel how many bits to use for each datum
     * @param minValue the minimum value of the data (used for normalization)
     * @param maxValue the maximum value of the data (used for normalization)
     */
    public PacketPacker(int bitsPerChannel, double minValue, double maxValue) {
        mBitDepth = bitsPerChannel;
        mMinValue = minValue;
        mRange = maxValue - minValue;
    }

    /**
     * Construct a new {@code PacketPacker}
     *
     * @param bitsPerChannel how many bits to use for each datum
     * @param maxValue the maximum value of the data (used for normalization)
     */
    public PacketPacker(int bitsPerChannel, double maxValue) {
        mBitDepth = bitsPerChannel;
        mMinValue = 0;
        mRange = maxValue;
    }

    /**
     * Pack a data array into a byte packet based on the parameters assigned to this
     * {@code PacketPacker}
     * <p>
     *     The data is normalized based on the min/max value supplied during construction
     * </p>
     *
     * @param data the data to normalize and pack
     *
     * @return the packed data (array of bytes)
     */
    public byte[] pack(float[] data) {
        int numChannels = data.length;
        BitSet bits = new BitSet(mBitDepth *numChannels);
        for(int i = 0; i < numChannels; i++) {
            BitSet singleChannelBits = toBitSet((data[i] - mMinValue)/mRange, mBitDepth);
            for(int b = 0; b < mBitDepth; b++) {
                int bitIndex = i* mBitDepth + b;
                bits.set(bitIndex, singleChannelBits.get(b));
            }
        }
        return Arrays.copyOf(toByteArray(bits), (int)Math.ceil(mBitDepth *numChannels/8.0));
    }

    /**
     * Pack a data array into a byte packet
     * <p>
     *     Assumes the data is already normalized
     * </p>
     *
     * @param data the normalized data to pack
     * @param bitDepth the number of bits used to pack each data point
     *
     * @return the packed data (array of bytes with length = data.length*bitDepth/8)
     */
    public static byte[] pack(float[] data, int bitDepth) {
        int numChannels = data.length;
        BitSet bits = new BitSet(bitDepth*numChannels);
        for(int i = 0; i < numChannels; i++) {
            BitSet singleChannelBits = toBitSet(data[i], bitDepth);
            for(int b = 0; b < bitDepth; b++) {
                int bitIndex = i*bitDepth + b;
                bits.set(bitIndex, singleChannelBits.get(b));
            }
        }
		return Arrays.copyOf(toByteArray(bits), (int)Math.ceil(bitDepth*numChannels/8.0));
    }

    // Convert a BitSet to an array of bytes
    private static byte[] toByteArray(BitSet bits) {
		byte[] bytes = new byte[(bits.length()+7)/8];
        for(int n = 0; n < 8 * bytes.length; n++) {
            int val = bits.get(n) ? 1 : 0;
            bytes[n/8] |= val<<(n%8);
        }
        return bytes;
    }

    // convert a float to a BitSet,
    private static BitSet toBitSet(double d, int bitDepth) {
        BitSet bits = new BitSet(bitDepth);
        long L = (long) (d * (Math.pow(2, bitDepth) - 1));
        int index = 0;
        while(L != 0) {
            if (L % 2L != 0) {
                bits.set(index);
            }
            ++index;
            L = L >>> 1;
        }
        return bits;
    }
}
