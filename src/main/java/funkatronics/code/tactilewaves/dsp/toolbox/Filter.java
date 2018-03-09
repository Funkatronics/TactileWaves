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
import java.util.logging.Logger;

/**
 * Class for filtering data series using recursive filters
 * <p>
 *     Butterworth, Critically Damped, and Bessel filter types are supported via instantiation with
 *     an arbitrary number of poles. A static method for filtering using any provided filter
 *     coefficients is available for custom filtering.
 * </p>
 *
 * @author Marco Martinez
 */

public class Filter {

    private static final Logger LOG = Logger.getLogger(Filter.class.getSimpleName());

    // PI constant, because it is used a lot in this class
    private final static double pi = Math.PI;

    // Integer constants that respresent the different filter types
    /**
     * A Butterworth filter
     */
    public static final int FILTER_BUTTERWORTH = 0;

    /**
     * A Critically Damped filter
     */
    public static final int FILTER_CRITICALLY_DAMPED = 1;

    /**
     * A Bessel filter
     */
    public static final int FILTER_BESSEL = 2;

    // Is this filter a high pass filter?
    private boolean mHighPass = false;

    // What filter type is this filter? (default = Butterworth)
    private int mFilterType = FILTER_BUTTERWORTH;

    // How many poles in this filter? (default = 2)
    private int mNP = 2;

    // What is the cutoff frequency of this filter? (scaled 0.0-0.5)
    private double mFc;

    // Arrays to hold filter coefficients
    private double[] a = new double[3];
    private double[] b = new double[3];

    /**
     * Construct a default filter
     *
     * @param Fc the cutoff frequency, as a fraction of the sample rate (0.0-0.5)
     * @param highPass is this filter a high pass filter?
     */
    public Filter(float Fc, boolean highPass){
        mFc = Fc;
        mNP = 2;
        mFilterType = FILTER_BUTTERWORTH;
        mHighPass = highPass;
        getCoefficients();
    }

    /**
     * Construct a new {@code Filter} object with these parameters
     *
     * @param Fc the cutoff frequency, as a fraction of the sample rate (0.0-0.5)
     * @param NP the number of poles of the filter
     * @param filterType the type of filter to use
     * @param highPass is this filter a high pass filter?
     */
    public Filter(float Fc, int NP, int filterType, boolean highPass){
        mFc = Fc;
        mNP = NP;
        mFilterType = filterType;
        mHighPass = highPass;
        getCoefficients();
    }

    // Get the filter coefficients for this filter instance
    private void getCoefficients() {
        if(mNP < 2) throw new IllegalArgumentException("Filter cannot have less than 2 poles!");
        double[] A = new double[mNP +3];
        double[] B = new double[mNP +3];
        A[2] = 1;
        B[2] = 1;
        for(int p = 1; p <= mNP /2; p++) {

            //gcHelper_old(Fc, NP, p, highPass);
            get2PoleCoefficients();

            double[] TA = Arrays.copyOf(A, A.length);
            double[] TB = Arrays.copyOf(B, B.length);

            for(int i = 2; i < A.length; i++) {
                A[i] = a[0]*TA[i] + a[1]*TA[i-1] + a[2]*TA[i-2];
                B[i] = TB[i] - b[1]*TB[i-1] - b[2]*TB[i-2];
            }
        }

        B[2] = 0;
        for(int i = 0; i < A.length-2; i++) {
            A[i] = A[i+2];
            B[i] = -B[i+2];
        }

        double SA = 0;
        double SB = 0;
        for(int i = 0; i < A.length-2; i++) {
            if(mHighPass) {
                SA += A[i]*Math.pow(-1, i);
                SB += B[i]*Math.pow(-1, i);
            } else {
                SA += A[i];
                SB += B[i];
            }
        }

        double gain = SA/(1 - SB);

        for(int i = 0; i < A.length-2; i++) {
            A[i] = A[i]/gain;
        }

        a = Arrays.copyOf(A, mNP +1);
        b = Arrays.copyOf(B, mNP +1);
    }

    // Helper function for getCoefficients() to get the 2-pole filter coefficients
    private void get2PoleCoefficients() {
        int numPass = 1;
        double c, g, p;
        switch (mFilterType) {
            case FILTER_BUTTERWORTH:
                c = Math.pow((Math.pow(2.0, 1.0/(double)numPass) - 1), -0.25);
                g = 1.0;
                p = Math.sqrt(2.0);
                break;
            case FILTER_CRITICALLY_DAMPED:
                c = Math.pow((Math.pow(2.0, 1.0/(2.0*(double)numPass)) - 1), -0.5);
                g = 1.0;
                p = 2.0;
                break;
            case FILTER_BESSEL:
                c = Math.pow((Math.pow((Math.pow(2.0, 1.0/(double)numPass) - 0.75), 0.5) - 0.5), -0.5) / Math.sqrt(3.0);
                g = 3.0;
                p = 3.0;
                break;
            default:
                c = Math.pow((Math.pow(2.0, 1 / numPass) - 1), -0.25);
                g = 1.0;
                p = Math.sqrt(2.0);
                mFilterType = FILTER_BUTTERWORTH;
        }

        if(mHighPass) c = 1/c;

        double fc = c * mFc;

        // Check stability
        if (!(fc > 0.0 && fc < 0.25))
            LOG.info("Filter unstable for Fc = " + mFc);
            //Log.i(TAG,"Filter unstable: " + fc);

        if(mHighPass) fc = 0.5 - fc;

        double w0 = Math.tan(pi*fc);

        double K1 = p*w0;
        double K2 = g*w0*w0;

        a[0] = K2/(1 + K1 + K2);
        a[1] = 2* a[0];
        a[2] = a[0];
        b[1] = 2* a[0]*(1/K2 - 1);
        b[2] = 1 - (a[0] + a[1] + a[2] + b[1]);

        if(mHighPass) {
            a[1] = -a[1];
            b[1] = -b[1];
        }
    }

    private void gcHelper_old(double Fc, int NP, int p, boolean highPass) {
        double rp = -Math.cos(pi/(NP*2) + (p-1)*pi/NP);
        double ip = Math.sin(pi/(NP*2) + (p-1)*pi/NP);

        double T = 2*Math.tan(0.5);
        double w = 2*pi*Fc;
        double M = rp*rp + ip*ip;
        double D = 4 - 4*rp*T + (M*T*T);
        double X0 = T*T/D;
        double X1 = 2*T*T/D;
        double X2 = T*T/D;
        double Y1 = (8 - 2*M*T*T)/D;
        double Y2 = (-4 - 4*rp*T - M*T*T)/D;

        double K = Math.sin(0.5 - w/2)/Math.sin(0.5 + (w/2));
        if(highPass) K = -Math.cos(w/2 + 0.5)/Math.cos(w/2 - 0.5);
        D = 1.0 + Y1*K - Y2*K*K;
        a[0] = (X0 - X1*K + X2*K*K)/D;
        a[1] = (-2*X0*K + X1 + X1*K*K - 2*X2*K)/D;
        a[2] = (X0*K*K - X1*K + X2)/D;
        b[1] = (2*K + Y1 + Y1*K*K - 2*Y2*K)/D;
        b[2] = (-(K*K) - Y1*K + Y2)/D;
        if(highPass) {
            a[1] = -a[1];
            b[1] = -b[1];
        }
    }

    private void getCoefficients_old() {
        int numPass = 2;
        double c, g, p;
        switch (mFilterType) {
            case FILTER_BUTTERWORTH:
                c = Math.pow((Math.pow(2.0, 1.0/(double)numPass) - 1), -0.25);
                g = 1.0;
                p = Math.sqrt(2.0);
                break;
            case FILTER_CRITICALLY_DAMPED:
                c = Math.pow((Math.pow(2.0, 1.0/(2.0*(double)numPass)) - 1), -0.5);
                g = 1.0;
                p = 2.0;
                break;
            case FILTER_BESSEL:
                c = Math.pow((Math.pow((Math.pow(2.0, 1.0/(double)numPass) - 0.75), 0.5) - 0.5), -0.5) / Math.sqrt(3.0);
                g = 3.0;
                p = 3.0;
                break;
            default:
                c = Math.pow((Math.pow(2.0, 1 / numPass) - 1), -0.25);
                g = 1.0;
                p = Math.sqrt(2.0);
                mFilterType = FILTER_BUTTERWORTH;
        }

        if(mHighPass) c = 1/c;

        double fc = c * mFc;

        // Check stability
        if (!(fc > 0.0 && fc < 0.25))
            System.out.println("Filter unstable: " + fc);

        if(mHighPass) fc = 0.5 - fc;

        double w0 = Math.tan(Math.PI*fc);

        double K1 = p*w0;
        double K2 = g*w0*w0;

        a[0] = K2/(1 + K1 + K2);
        a[1] = 2* a[0];
        a[2] = a[0];
        b[1] = 2* a[0]*(1/K2 - 1);
        b[2] = 1 - (a[0] + a[1] + a[2] + b[1]);

        if(mHighPass) {
            a[1] = -a[1];
            b[1] = -b[1];
        }
    }

    /**
     * Enable or disable the high pass flag
     *
     * @param highPass true to enable high pass filtering, false for low pass
     */
    public void enableHighPass(boolean highPass) {
        mHighPass = highPass;
        getCoefficients();
    }

    /**
     * Sets the filter type of this fiter
     *
     * @param filterType the type of filter to use, chosen from {@code FILTER_BUTTERWORTH},
     * {@code FILTER_CRITICALLY_DAMPED}, or {@code FILTER_BESSEL}.
     */
    public void setFilterType(int filterType) {
        mFilterType = filterType;
        getCoefficients();
    }

    /**
     * Filter a signal/data series
     *
     * @param x the signal or data series to filtered
     *
     * @return the filtered data
     */
    public float[] filter(float[] x) {
        return filter(a, b, x);
    }

    /**
     * Filter a signal/data buffer using the supplied filter coefficients
     *
     * @param a the numerator filter coefficients
     * @param b the denominator filter coefficients
     * @param x the signal or data series to be filtered
     *
     * @return the filtered data
     */
    public static float[] filter(float[] a, float[] b, float[] x) {
        if(a.length < 1 || b.length < 1)
            throw new IllegalArgumentException("Filter coefficients cannot be empty");
        float[] y = new float[x.length];
        for(int n = 0; n < x.length; n++) {
            for (int k = 0; k < a.length; k++)
                if(k <= n) y[n] += a[k] * x[n - k];
            for (int k = 1; k < b.length; k++)
                if(k <= n) y[n] += b[k] * y[n - k];
        }
        return y;
    }

    /**
     * Filter a signal/data buffer using the supplied filter coefficients
     *
     * @param a the numerator filter coefficients
     * @param b the denominator filter coefficients
     * @param x the signal or data series to be filtered
     *
     * @return the filtered data
     */
    public static float[] filter(double[] a, double[] b, float[] x) {
        if(a.length < 1 || b.length < 1)
            throw new IllegalArgumentException("Filter coefficients cannot be empty");
        float[] y = new float[x.length];
        for(int n = 0; n < x.length; n++) {
            for (int k = 0; k < a.length; k++)
                if(k <= n) y[n] += a[k] * x[n - k];
            for (int k = 1; k < b.length; k++)
                if(k <= n) y[n] += b[k] * y[n - k];
        }
        return y;
    }

    /**
     * Filter a signal/data buffer using the supplied filter coefficients
     *
     * @param a the numerator filter coefficients
     * @param b the denominator filter coefficients
     * @param x the signal or data series to be filtered
     *
     * @return the filtered data
     */
    public static double[] filter(double[] a, double[] b, double[] x) {
        if(a.length < 1 || b.length < 1)
            throw new IllegalArgumentException("Filter coefficients cannot be empty");
        double[] y = new double[x.length];
        for(int n = 0; n < x.length; n++) {
            for (int k = 0; k < a.length; k++)
                if(k <= n) y[n] += a[k] * x[n - k];
            for (int k = 1; k < b.length; k++)
                if(k <= n) y[n] += b[k] * y[n - k];
        }
        return y;
    }

    /**
     * Apply a Pre-Emphasis Filter toa signal/data series
     *
     * @param x the signal or data series to be filtered
     *
     * @return the filtered data
     */
    public static float[] preEmphasis(float[] x) {
        double[] a = {1.0, -0.63};
        double[] b = {1.0};
        return Filter.filter(b, a, x);
    }

    /**
     * Apply a Pre-Emphasis Filter toa signal/data series
     *
     * @param x the signal or data series to be filtered
     *
     * @return the filtered data
     */
    public static double[] preEmphasis(double[] x) {
        double[] a = {1.0, -0.63};
        double[] b = {1.0};
        return Filter.filter(b, a, x);
    }

    /**
     * Apply a Moving Average Filter to a signal/data series
     *
     * @param x the signal or data series to be filtered
     * @param windowSize the size of the moving average window
     *
     * @return the filtered data
     */
    public static float[] movingAvg(float[] x, int windowSize) {
        double[] b = new double[windowSize];
        for(int i = 0; i < windowSize; i++)
            b[i] = 1.0/windowSize;
        double[] a = new double[1];
        a[0] = 1.0;
        return Filter.filter(b, a, x);
    }

    /**
     * Apply a Moving Average Filter to a signal/data series
     *
     * @param x the signal or data series to be filtered
     * @param windowSize the size of the moving average window
     *
     * @return the filtered data
     */
    public static double[] movingAvg(double[] x, int windowSize) {
        double[] b = new double[windowSize];
        for(int i = 0; i < windowSize; i++)
            b[i] = 1.0/windowSize;
        double[] a = new double[1];
        a[0] = 1.0;
        return Filter.filter(b, a, x);
    }

    private static void LPWindowedSincKernel(float[] h, float fc) {
        int M = h.length - 1;
        float sum = 0;
        for(int i = 0; i <= M; i++) {
            if(i == M/2) h[i] = (float)(2*pi*fc);
            else h[i] = (float)(Math.sin(2*pi*fc*(i - (M/2)))/(i - (M/2)));
            h[i] *= (0.54 - 0.46*Math.cos(2*pi*i/M));
            sum += h[i];
        }
        for(int i = 0; i <= M; i++) {
            h[i] /= sum;
        }
    }

    private static void LPWindowedSincKernel(double[] h, float fc) {
        int M = h.length - 1;
        double sum = 0;
        for(int i = 0; i <= M; i++) {
            if(i == M/2) h[i] = 2*pi*fc;
            else h[i] = Math.sin(2*pi*fc*(i - (M/2)))/(i - (M/2));
            h[i] *= 0.54 - 0.46*Math.cos(2*pi*i/M);
            sum += h[i];
        }
        for(int i = 0; i <= M; i++) {
            h[i] /= sum;
        }
    }
}
