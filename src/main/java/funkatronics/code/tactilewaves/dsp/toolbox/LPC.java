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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import funkatronics.code.tactilewaves.dsp.utilities.Complex;
import funkatronics.code.tactilewaves.dsp.utilities.RootSolver;
import funkatronics.code.tactilewaves.dsp.utilities.SolverNotConvergedException;


/**
 * Class for computing the forward Linear Prediction Coefficients (LPC) of a data series.
 * <p>
 *     This class is used primarily to estimate the formant frequencies of spoken audio, but a
 *     generalized {@code LPC} method is provided for general LPC usage.
 * </p>
 * <p>
 *     Levinson-Durbin Recursion is used to solve the Youle-Walker equations and compute the linear
 *     predictor coefficients.
 * </p>
 *
 * @author Marco Martinez
 */

public class LPC {

    // do not allow instantiation
    private LPC() {}

    /**
     * Estimate the formant frequencies of an audio signal
     *
     * @param x the audio signal
     * @param numFormants the number of formants to estimate
     * @param Fs the sampleRate of the audio signal
     * @return the formant frequencies and bandwidths in a 2D array
     *         ({@code {{Formant 1, Bandwidth 1}, {Formant 2, Bandwidth 2}, etc.}})
     */
    public static double[][] estimateFormants(float[] x, int numFormants, int Fs) {
        double[] result = LPC.LPC(x, 2*numFormants + 2);
        try {
            Complex[] p = RootSolver.roots(result, false);
            return LPC.LPCRoots2Formants(p,Fs);
        } catch (SolverNotConvergedException e) {
            return new double[][]{};
        }
    }

    /**
     * Estimate the formant frequencies of an audio signal
     *
     * @param x the audio signal
     * @param numFormants the number of formants to estimate
     * @param Fs the sampleRate of the audio signal
     * @return the formant frequencies and bandwidths in a 2D array
     *         ({@code {{Formant 1, Bandwidth 1}, {Formant 2, Bandwidth 2}, etc.}})
     */
    public static double[][] estimateFormants(double[] x, int numFormants, int Fs) {
        double[] result = LPC.LPC(x, 2*numFormants + 2);
        try {
            Complex[] p = RootSolver.roots(result, false);
            return LPC.LPCRoots2Formants(p,Fs);
        } catch (SolverNotConvergedException e) {
            return new double[][]{};
        }
    }

//    public static float[] autocorr(float[] x, int lags) {
//        float[] ac = new float[lags];
//        float max = -Float.MAX_VALUE;
//        for(int tau = 0; tau < lags; tau++) {
//            ac[tau] = autocorrAtTau(x, tau);
//            if(ac[tau] > max) max = ac[tau];
//        }
//        for(int i = 0; i < ac.length; i++) {
//            ac[i] /= max;
//        }
//        return ac;
//    }
//
//    public static double[] autocorr(double[] x, int lags) {
//        double[] ac = new double[lags];
//        double max = -Double.MAX_VALUE;
//        for(int tau = 0; tau < lags; tau++) {
//            ac[tau] = autocorrAtTau(x, tau);
//            if(ac[tau] > max) max = ac[tau];
//        }
//        for(int i = 0; i < ac.length; i++) {
//            ac[i] /= max;
//        }
//        return ac;
//    }
//
//    public static float[] autocorrFFT(float[] x) {
//        return FFT.autocorr(x);
//    }
//
//    public static double[] autocorrFFT(double[] x) {
//        return FFT.autocorr(x);
//    }
//
//    public static float autocorrAtTau(float[] x, int tau) throws IndexOutOfBoundsException{
//        float ac = 0.0f;
//        for(int j = tau; j < x.length; j++){
//            ac += x[j] * x[j-tau];
//        }
//        return ac;
//    }
//
//    public static double autocorrAtTau(double[] x, int tau) throws IndexOutOfBoundsException{
//        double ac = 0.0f;
//        for(int j = tau; j < x.length; j++){
//            ac += x[j] * x[j-tau];
//        }
//
//        return ac;
//    }

    /**
     * Get coefficients of a pth-order linear predictor (FIR Filter) of a data series
     *
     * @param x the data series to poly-fit
     * @param order the oder of the prediction filter
     * @return the ordered coefficients of the linear predictor
     */
    public static double[] LPC(float[] x, int order) throws IllegalArgumentException {
        if(order > x.length)
            throw new IllegalArgumentException("x must be an array with length greater or equal to the prediction order");
        double[] output = new double[order];

        double[] xd = new double[x.length];
        for(int i = 0; i < x.length; i++) xd[i] = (double)x[i];

        double[] R = FFT.autocorr(xd);

        output = levinson(R, order, true);

        return output;
    }

    /**
     * Get coefficients of a pth-order linear predictor (FIR Filter) of a data series
     *
     * @param x the data series to poly-fit
     * @param order the oder of the prediction filter
     * @return the ordered coefficients of the linear predictor
     */
    public static double[] LPC(double[] x, int order) throws IllegalArgumentException {
        if(order > x.length)
            throw new IllegalArgumentException("x must be an array with length greater or equal to the prediction order");
        double[] output = new double[order];

        double[] R = FFT.autocorr(x);

        output = levinson(R, order, true);

        return output;
    }

    // Implementation of Levinson-Durbin Recursion
    private static double[] levinson(double[] r, int order, boolean allow_singularity) throws IllegalArgumentException {
        double T0 = r[0];
        double[] T = Arrays.copyOfRange(r, 1, r.length);
        int M;

        if (order == -1)
            M = T.length;
        else
            M = order;

        double[] A = new double[M];
        double[] ref = new double[M];

        double P = T0;

        for (int k = 0; k < M; k++) {
            double save = T[k];
            double temp;

            if (k == 0) temp = -save / P;
            else {
                for (int j = 0; j < k; j++) save += A[j] * T[k - j - 1];
                temp = -save / P;
            }
            P = P * (1. - Math.pow(temp, 2));
            if (P <= 0 && !allow_singularity)
                throw new IllegalArgumentException("Singular matrix, but singularity is not allowed");

            A[k] = temp;
            // Save reflection coefficient
            ref[k] = temp;
            if (k == 0)
                continue;

            int khalf = (k + 1) / 2;
            for (int j = 0; j < khalf; j++) {
                int kj = k - j - 1;
                save = A[j];
                A[j] = save + temp * A[kj];
                if (j != kj)
                    A[kj] += temp * save;
            }

        }
        double[] result = new double[A.length+1];
        result[0] = 1.0;
        System.arraycopy(A, 0, result, 1, A.length);
        return result;
    }

    /**
     * Convert the roots of a polynomial to formant frequencies
     *
     * @param r the array of complex polynomial roots
     * @param Fs the sample rate to convert frequencies
     * @return the formant frequencies and bandwidths in a 2D array
     *         ({@code {{Formant 1, Bandwidth 1}, {Formant 2, Bandwidth 2}, etc.}})
     */
    public static double[][] LPCRoots2Formants(Complex[] r, int Fs) {
        ArrayList<double[]> result = new ArrayList<>();
        double fmax = Fs/2;
        for (Complex c : r) {
            double real = c.real();
            double imag = c.imag();
            if (imag >= 0) {
                double freq = Math.atan2(imag, real) * ((double) Fs / (2 * Math.PI));
                double bw = (-0.5) * (Fs / (2 * Math.PI)) * Math.log(c.mag());
                if (freq < fmax) result.add(new double[]{freq, bw});
            }
        }
        if(result.isEmpty()) return new double[][] {{0.0,0.0}};

        // Sorting with Lambda function (Java 8+)
//        result.sort((double[] d1, double[] d2)->Double.compare(d1[0], d2[0]));

        double[][] sorted = result.toArray(new double[result.size()][2]);

        // Sorting with Lambda function (Java 8+)
//        Arrays.sort(sorted, (double[] d1, double[] d2)->Double.compare(d1[0], d2[0]));

        Arrays.sort(sorted,new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                if(o1[0] < o2[0]) return -1;
                else if(o1[0] > o2[0]) return 1;
                else return 0;
            }
        });
        return sorted;
//        return result.toArray(new double[2][result.size()]);
    }
}
