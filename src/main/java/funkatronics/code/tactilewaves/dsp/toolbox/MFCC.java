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

package funkatronics.code.tactilewaves.dsp.toolbox;

import java.util.Arrays;

/**
 * Class for computing the Mel Frequency Cepstrum Coefficients of an audio signal.
 *
 * @author Marco Martinez
 */

public class MFCC {

    // INteger constants that define the type of litering to use
    public final static int LIFTERING_NONE = -1;
    public final static int LIFTERING_LINEAR = 0;
    public final static int LIFTERING_SINUSOIDAL = 1;
    public final static int LIFTERING_EXPONENTIAL = 2;

    // Do not allow instantiation
    private MFCC() {}

    /**
     * Get the Mel Frequency Cepstrum Coefficients of an audio signal
     * <p>
     *     This uses 26 filter banks from 0 Hz to Fs/2 and returns 13 MFCCs
     * </p>
     *
     * @param x the audio signal
     * @param Fs the sample rate of the audio signal
     *
     * @return the Mel Frequency Cepstrum Coefficients
     */
    public static double[] getMFCCs(float[] x, int Fs) {
        return getMFCCs(x, Fs, 0, Fs/2, 26, 13, LIFTERING_SINUSOIDAL);
    }

    /**
     * Get the Mel Frequency Cepstrum Coefficients of an audio signal
     * <p>
     *     This uses 26 filter banks from 0 Hz to Fs/2 and returns 13 MFCCs
     * </p>
     *
     * @param x the audio signal
     * @param Fs the sample rate of the audio signal
     *
     * @return the Mel Frequency Cepstrum Coefficients
     */
    public static double[] getMFCCs(double[] x, int Fs) {
        return getMFCCs(x, Fs, 0, Fs/2, 26, 13, LIFTERING_SINUSOIDAL);
    }

    /**
     * Get the Mel Frequency Cepstrum Coefficients of an audio signal
     *
     * @param x the audio signal
     * @param Fs the sample rate of the audio signal
     * @param minF the minimum frequency of the filter bank
     * @param maxF the maximum frequency of the filter bank
     * @param Nfilt the number of bands in the filter bank
     * @param Ncoeffs the number of cepstrum coefficients to return
     * @param lifter the type of liftering to use on the computed coefficients. Choose from
     *               {@code LIFTERING_NONE}, {@code LIFTERING_LINEAR}, {@code LIFTERING_SINUSOIDAL},
     *               and {@code LIFTERING_EXPONENTIAL}.
     *
     * @return the Mel Frequency Cepstrum Coefficients
     */
    public static double[] getMFCCs(float[] x, int Fs, float minF, float maxF, int Nfilt,
                                    int Ncoeffs, int lifter) {
        if(minF < 0 || minF > maxF)
            throw new IllegalArgumentException("Minimum frequency out of bounds: 0 <= minF < maxF");
        if(maxF > Fs/2)
            throw new IllegalArgumentException("Maximum frequency out of bounds: minF < maxF <= Fs/2");
        if(Ncoeffs > Nfilt)
            throw new IllegalArgumentException("Cannot return more coefficients than there are filter bands: Ncoeffs <= Nfilt");

        // Get Periodogram
        float[] X = FFT.PowerSpectrum(x);

        // Initialize filter bank
        int[] f = initFilterBanks(Nfilt, minF, maxF, X.length, Fs);

        // Compute filter bank
        double[] fb = new double[Nfilt];
        for(int m = 1; m <= Nfilt; m++) {
            for(int k = f[m-1]; k <= f[m]; k++) fb[m-1] += X[k] * (k - f[m-1])/(f[m] - f[m-1]);
            for(int k = f[m]; k <= f[m+1]; k++) fb[m-1] += X[k] * (f[m+1] - k)/(f[m+1] - f[m]);
            //fb[m-1] = 20.0 * Math.log10(fb[m-1]);
            if(fb[m-1] < 1e-5 || fb[m-1] != fb[m-1]) fb[m-1] = 1e-5;
            fb[m-1] = Math.log(fb[m-1]);
        }

        // DCT to decorrelate (whiten) the filter coefficients
        DCT.DCT(fb, true);
        double[] MFCCs = Arrays.copyOfRange(fb, 1, Ncoeffs+1);

        switch(lifter) {
            case LIFTERING_LINEAR:
                linLifter(MFCCs);
                break;
            case LIFTERING_SINUSOIDAL:
                sinLifter(MFCCs);
                break;
            case LIFTERING_EXPONENTIAL:
                expLifter(MFCCs);
                break;
            default:
                // no liftering
                break;
        }

        return MFCCs;
    }

    /**
     * Get the Mel Frequency Cepstrum Coefficients of an audio signal
     *
     * @param x the audio signal
     * @param Fs the sample rate of the audio signal
     * @param minF the minimum frequency of the filter bank
     * @param maxF the maximum frequency of the filter bank
     * @param Nfilt the number of bands in the filter bank
     * @param Ncoeffs the number of cepstrum coefficients to return
     * @param lifter the type of liftering to use on the computed coefficients. Choose from
     *               {@code LIFTERING_NONE}, {@code LIFTERING_LINEAR}, {@code LIFTERING_SINUSOIDAL},
     *               and {@code LIFTERING_EXPONENTIAL}.
     *
     * @return the Mel Frequency Cepstrum Coefficients
     */
    public static double[] getMFCCs(double[] x, int Fs, float minF, float maxF, int Nfilt,
                                    int Ncoeffs, int lifter) {
        if(minF < 0 || minF > maxF)
            throw new IllegalArgumentException("Minimum frequency out of bounds: 0 <= minF < maxF");
        if(maxF > Fs/2)
            throw new IllegalArgumentException("Maximum frequency out of bounds: minF < maxF <= Fs/2");
        if(Ncoeffs > Nfilt)
            throw new IllegalArgumentException("Cannot return more coefficients than there are filter bands: Ncoeffs <= Nfilt");

        // Get Periodogram
        double[] X = FFT.PowerSpectrum(x);

        // Initialize filters
        int[] f = initFilterBanks(Nfilt, minF, maxF, X.length, Fs);

        // Compute filter bank
        double[] fb = new double[Nfilt];
        for(int m = 1; m <= Nfilt; m++) {
            for(int k = f[m-1]; k <= f[m]; k++) fb[m-1] += X[k] * (k - f[m-1])/(f[m] - f[m-1]);
            for(int k = f[m]; k <= f[m+1]; k++) fb[m-1] += X[k] * (f[m+1] - k)/(f[m+1] - f[m]);
            //fb[m-1] = 20.0 * Math.log10(fb[m-1]);
            if(fb[m-1] < 1e-5 || fb[m-1] != fb[m-1]) fb[m-1] = 1e-5;
            fb[m-1] = Math.log(fb[m-1]);
        }

        // DCT to decorrelate (whiten) the filter coefficients
        DCT.DCT(fb, true);
        double[] MFCCs = Arrays.copyOfRange(fb, 1, Ncoeffs+1);

        switch(lifter) {
            case LIFTERING_LINEAR:
                linLifter(MFCCs);
                break;
            case LIFTERING_SINUSOIDAL:
                sinLifter(MFCCs);
                break;
            case LIFTERING_EXPONENTIAL:
                expLifter(MFCCs);
                break;
            default:
                // no liftering
                break;
        }

        return MFCCs;
    }

    /**
     * Compute delta (differential) between two arrays
     * <p>
     *     Used with frames of MFCCs to compute the change between frames (differential
     *     coefficients).
     * </p>
     * @param T the current frame
     * @param TM1 the previous frame
     * @return the element wise delta (difference) between these two frames
     */
    public static double[] deltas(double[] T, double[] TM1) {
        if(T.length != TM1.length)
            throw new IllegalArgumentException("Current frame and previous frame must be the same length");
        double[] deltas = new double[T.length];
        for(int i = 0; i < deltas.length; i++)
            deltas[i] = (T[i] - TM1[i]);
        return deltas;
    }

    /**
     * Initialize the Mel Filter Bank
     * @param Nfilt the number of individual filter bands in the filter bank
     * @param minF the minimum frequency (Hz) of the filter bank
     * @param maxF the minimum frequency (Hz) of the filter bank
     * @param Nfft the size of the incoming audio/FFT buffer
     * @param Fs the sample rate of the incoming audio
     *
     * @return an array of indices representing the filter bank
     */
    public static int[] initFilterBanks(int Nfilt, float minF, float maxF, int Nfft, int Fs) {
        double[] m = new double[Nfilt + 2];
        m[0] = freq2Mel(minF);
        m[m.length-1] = freq2Mel(maxF);

        for(int i = 1; i <= Nfilt; i++) m[i] = m[0] + i*(m[m.length-1] - m[0])/(Nfilt+1);

        for(int i = 0; i < m.length; i++) m[i] = mel2Freq(m[i]);

        int[] f = new int[Nfilt + 2];
        for(int i = 0; i < m.length; i++) f[i] = (int)Math.floor((Nfft+1)*m[i]/Fs);

        return f;
    }

    // Perform linear liftering on an array
    private static void linLifter(double[] array) {
        for(int i = 1; i <= array.length; i++)
            array[i-1] *= i;
    }

    // Perform sinusoidal liftering on an array
    private static void sinLifter(double[] array) {
        double D = array.length;
        for(int i = 1; i <= array.length; i++)
            array[i-1] *= (1.0 + (D/2.0)*Math.sin(i*Math.PI/D));
    }

    // Perform exponential liftering on an array
    private static void expLifter(double[] array) {
        double s = 1.5;
        double tau = 5.0;
        for(int i = 1; i <= array.length; i++)
            array[i-1] *= (Math.pow(i, s)*Math.exp(-i*i/(2*tau*tau)));
    }

    // Convert a frequency in Hz to mel scale
    private static float freq2Mel(float f) {
        return (float)freq2Mel((double)f);
    }

    // Convert a frequency in Hz to mel scale
    private static double freq2Mel(double f) {
        // maybe should use 1127 instead of 1125
        //return (1125.0 * Math.log(1 + f/700));
        return 2595.0 * Math.log10(1 + f/700.0);
    }

    // Convert a mel scale value to frequency (Hz)
    private static double mel2Freq(double mel) {
        // maybe should use 1127 instead of 1125
        //return 700.0 * (Math.exp(mel/1125.0) - 1);
        return 700.0 * (Math.pow(10, mel/2595.0) - 1);
    }

}
