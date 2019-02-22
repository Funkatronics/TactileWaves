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
     * Compute the real Cepstrum of a signal in-place
     *
     * @param x the signal to process
     */
    public static void rCepstrum(float[] x) {
        // Real FFT
        float[] ImX = FFT.fft(x);

        for(int i = 0; i < x.length; i++) {
            x[i] = (float) Math.log(Math.sqrt(x[i]*x[i] + ImX[i]*ImX[i]));
            ImX[i] = 0;
        }

        FFT.ifft(x, ImX);
    }

    /**
     * Compute the power Cepstrum of a signal in-place
     *
     * @param x the signal to process
     */
    public static void pCepstrum(float[] x) {
        // Real FFT
        float[] ImX = FFT.fft(x);

        for(int i = 0; i < x.length; i++) {
            x[i] = (float) Math.log(x[i]*x[i] + ImX[i]*ImX[i]);
            ImX[i] = 0;
        }

        FFT.ifft(x, ImX);

		for(int i = 0; i < x.length; i++) {
			x[i] = x[i]*x[i];
		}
    }

    /**
     * Compute the complex cepstrum of a signal in-place
     *
     * @param x the signal to process
	 *
	 * @return the imaginary part of the cepstrum
     */
    public static int cCepstrum(float[] x) {
        float[] ImX = FFT.fft(x);

		for(int i = 0; i < x.length; i++) {
            float temp = x[i];
            x[i] = (float) Math.log(Math.sqrt(x[i]*x[i] + ImX[i]*ImX[i]));
            ImX[i] = (float) Math.atan2(ImX[i],temp);
            // unwrap phase
            if(i > 0) {
                if((ImX[i] - ImX[i-1]) >= Math.PI) {
                    ImX[i] -= 2 * Math.PI;
                } else if((ImX[i] - ImX[i-1]) <= -Math.PI) {
                    ImX[i] += 2 * Math.PI;
                }
            }
        }
        int nd = rcUnwrap(ImX);

        FFT.ifft(x, ImX);

        return nd;
    }

	/**
	 * Compute the inverse complex cepstrum of a signal in-place
	 *
	 * @param x the real part of the cepstrum
	 */
	public static void icCepstrum(float[] x, int nd) {
		float[] ImX = FFT.fft(x);

		rcwrap(ImX, nd);
		for(int i = 0; i < x.length; i++) {
			// wrap phase
			if(i > 0) {
				if((ImX[i] - ImX[i-1]) >= Math.PI) {
					ImX[i] += 2 * Math.PI;
				} else if((ImX[i] - ImX[i-1]) <= -Math.PI) {
					ImX[i] -= 2 * Math.PI;
				}
			}
			double temp = Math.exp(x[i]);
			x[i] = (float) (temp*Math.cos(ImX[i]));
			ImX[i] = (float) (temp*Math.sin(ImX[i]));
		}

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
        float[] xc = Window.hamming(x);
        rCepstrum(xc);
        lifterHT(xc, Fs/500);
        int loc = Utilities.findHighestPeaks(xc, 1)[0];
        if(xc[loc] > 0.01f) return (float)Fs/loc;
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
		float[] xc = Window.hamming(x);
        rCepstrum(xc);
        lifterLT(xc, Fs/500);
        FFT.fft(xc);
        int[] locs = Utilities.findOrderedPeaks(xc, FNum);
        float[] formants = new float[FNum];
        for(int i = 0; i < FNum; i++) formants[i] = locs[i]*((float)Fs/xc.length);
        return formants;
    }

    /**
     * Perform Low-Time Liftering on an audio signal
     * <p>"Liftering" in the quefrency domain is equivalent to filtering in the frequency domain.
     * "Low-Time Liftering" is similar to low-pass filtering in that slowly varying components of
     * the signal are kept while quickly varying content is attenuated</p>
     *
     * @param c the Cepstrum to lifter
     * @param Lc the cuttoff length of the liftering window
     */
    public static void lifterLT(float[] c, int Lc) {
        for(int i = Lc; i < c.length; i++) {
            c[i] = 0;
        }
    }

    /**
     * Perform High-Time Liftering on an audio signal
     * <p>"Liftering" in the quefrency domain is equivalent to filtering in the frequency domain.
     * "High-Time Liftering" is similar to high-pass filtering in that quickly varying components of
     * the signal are kept while slowly varying content is attenuated</p>
     *
     * @param c the Cepstrum to lifter
     * @param Lc the cuttoff length of the liftering window
     */
    public static void lifterHT(float[] c, int Lc) {
        for(int i = 0; i < Lc; i++)
            c[i] = 0;
		for(int i = c.length/2; i < c.length; i++)
			c[i] = 0;
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
     */
    private static int rcUnwrap(float[] X) {
        int n = X.length;
        int nh = (n+1)/2;
        int idx = nh;
        if(n== 1) idx = 0;
        int nd = (int)Math.round(X[idx]/Math.PI);
        for(int i = 0; i < n; i++) {
            X[i] -= Math.PI*nd*i/nh;
        }
        return nd;
    }

	/**
	 * A special version of wrap that adds a straight line to the phase. This method is
	 * functionally identical to MatLab's rcwrap(), a local function from MatLab's own cepstral
	 * analysis package, <a href="https://www.mathworks.com/help/signal/ref/icceps.html">icceps</a>.
	 */
	private static void rcwrap(float[] X, int nd) {
		int n = X.length;
		int nh = (n+1)/2;
		for(int i = 0; i < n; i++) {
			X[i] += Math.PI*nd*i/nh;
		}
	}
}
