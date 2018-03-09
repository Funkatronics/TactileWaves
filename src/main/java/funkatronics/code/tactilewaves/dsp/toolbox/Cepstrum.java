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

import funkatronics.code.tactilewaves.dsp.utilities.Utilities;


/**
 * Class for performing Cepstral Transforms on audio data
 * <p>
 *     See <a href="http://iitg.vlab.co.in/?sub=59&brch=164&sim=615&cnt=1">this</a> website for a
 *     basic explanation of Cepstral Analysis.
 * </p>
 *
 * @author Marco Martinez
 * @see FFT
 * @see MFCC
 */

public class Cepstrum {

    // Do not allow instantiation
    private Cepstrum() {}

    /**
     * Compute the Real Cepstrum of a signal in-place
     *
     * @param x the signal to process
     */
    public static void rCepstrum(float[] x) {
        float[] ImX = new float[x.length];

        FFT.fft(x, ImX);

        for(int i = 0; i < x.length; i++) {
            x[i] = (float) Math.log(Math.sqrt(x[i]*x[i] + ImX[i]*ImX[i]));
            ImX[i] = 0;
        }

        FFT.ifft(x, ImX);
    }

    /**
     * Compute the Complex Cepstrum of a signal in-place
     *
     * @param x the signal to process
     */
    public static void cCepstrum(float[] x) {
        float[] ImX = new float[x.length];
        Arrays.fill(ImX, 0.0f);

        FFT.fft(x, ImX);

        for(int i = 0; i < x.length; i++) {
            float temp = x[i];
            x[i] = (float) (Math.log(Math.sqrt(x[i]*x[i] + ImX[i]*ImX[i])));
            ImX[i] = (float) (Math.atan2(ImX[i],temp));
            // unwrap phase
            if(i > 0) {
                if((ImX[i] - ImX[i-1]) >= Math.PI) {
                    ImX[i] -= 2 * Math.PI;
                } else if((ImX[i] - ImX[i-1]) <= -Math.PI) {
                    ImX[i] += 2 * Math.PI;
                }
            }
        }
        rcUnwrap(ImX);

        FFT.ifft(x, ImX);
    }

    /**
     * Estimate the pitch of an audio signal using Cepstral analysis and return it
     *
     * @param x the signal to process
     * @param Fs the sample rate of the audio (needed to convert pitch to Hz
     *
     * @return the estimated pitch in Hz
     */
    public static float estimatePitch(float[] x, int Fs) {
        Window.Hanning(x);
        rCepstrum(x);
        lifterLT(x, x.length/2);
        lifterHT(x, 15);
        int loc = Utilities.maxLoc(x);
        if(x[loc] > 0.07f) return (float)Fs/loc;
        else return -1.0f;
    }

    /**
     * Estimate the formant frequencies of an audio signal using Cepstral analysis
     *
     * @param x the signal to process
     * @param FNum the number of Formants to estimate
     * @param Fs the sample rate of the audio - used to convert formants to Hz.
     *
     * @return an array of estimated formant frequencies
     */
    public static float[] estimateFormants(float[] x, int FNum, int Fs) {
        Window.Hanning(x);
        rCepstrum(x);
        lifterLT(x, 15);
        float[] ImX = new float[x.length];
        FFT.fft(x, ImX);
        int[] locs = Utilities.findOrderedPeaks(x, FNum);
        float[] formants = new float[FNum];
        for(int i = 0; i < FNum; i++) formants[i] = locs[i]*Fs/x.length;
        return formants;
    }

    /**
     * Perform Low-Time Liftering on an audio signal
     * <p>"Liftering" in the quefrency domain is equivalent to filtering in the frequency domain.
     * "Low-Time Liftering" is similar to low-pass filtering in that slowly varying components of
     * the signal are kept while quickly varying content is attenuated</p>
     *
     * @param X the Cepstrum to lifter
     * @param Lc the cuttoff length of the liftering window
     */
    public static void lifterLT(float[] X, int Lc) {
        for(int i = Lc; i < X.length; i++) {
            X[i] = 0;
        }
    }

    /**
     * Perform High-Time Liftering on an audio signal
     * <p>"Liftering" in the quefrency domain is equivalent to filtering in the frequency domain.
     * "High-Time Liftering" is similar to high-pass filtering in that quickly varying components of
     * the signal are kept while slowly varying content is attenuated</p>
     *
     * @param X the Cepstrum to lifter
     * @param Lc the cuttoff length of the liftering window
     */
    public static void lifterHT(float[] X, int Lc) {
        for(int i = 0; i < Lc; i++) {
            X[i] = 0;
        }
    }

    /**
     * Unwraps a signal's phase by adding/subtracting appropriate multiples of 2π form its elements.
     *
     * @param X the signal to unwrap
     */
    private static void unwrap(float[] X) {
        for(int i = 0; i < X.length; i++) {
            if(i > 0) {
                if((X[i] - X[i-1]) >= Math.PI) {
                    X[i] -= 2 * Math.PI;
                } else if((X[i] - X[i-1]) <= -Math.PI) {
                    X[i] += 2 * Math.PI;
                }
            }
        }
    }

    /**
     * A special version of unwrap that subtracts a straight line from the phase. This method is
     * functionally identical to MatLab's rcunwrap(), a local function from MatLab's own cepstral
     * analysis package, <a href="https://www.mathworks.com/help/signal/ref/cceps.html">cceps</a>.
     *
     * @param X the signal to unwrap
     */
    private static void rcUnwrap(float[] X) {
        int n = X.length;
        unwrap(X); // unwrap the phase
        int nh = (n+1)/2;
        int idx = nh;
        if(X.length == 1) idx = 0;
        double nd = (X[idx]/Math.PI);
        for(int i = 0; i < X.length; i++) {
            X[i] = X[i] - (float) (Math.PI*nd*i/nh);
        }
    }
}
