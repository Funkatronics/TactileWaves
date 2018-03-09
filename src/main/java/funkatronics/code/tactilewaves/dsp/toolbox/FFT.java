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

import funkatronics.code.tactilewaves.dsp.utilities.Complex;

/**
 *  Class for performing Fast Fourier Transforms on sample buffers
 *  <p>
 *      Both {@link Complex} objects and single/double precision arrays are supported for radix-2 FFT's,
 *      as well as lookup based fft/ifft for performance (requires instantiation)
 *
 *      Any single/double precision array length is supported for Fourier transformation.
 *      A radix-2 fft will always be used if possible, applying Bluestein's algorithm otherwise.
 *  </p>
 *
 * @author Marco Martinez
 * @see Complex
 * @see Cepstrum
 * @see DCT
 * @see LPC
 * @see MFCC
 */
public class FFT {

    // The FFT length
    private int mSize;

    // Lookup tables. Only need to recompute when fft size is changed.
    private double[] mCos;
    private double[] mSin;

    /**
     * Construct a new {@code FFT} instance with the specified size
     *
     * @param size the length/size of the FFTs performed by this {@code FFT} instance
     */
    public FFT(int size) {
        // Make sure n is a power of 2
        if (size == 0 || (size&(size-1)) != 0)
            throw new RuntimeException("Length size must be power of 2");
        this.mSize = size;
        int M = (int)(Math.log(size)/Math.log(2));

        // precompute tables
        this.mCos = new double[M];
        this.mSin = new double[M];

        int N1 = 1;
        for(int i = 0; i < M; i++) {
            this.mCos[i] = Math.cos(Math.PI / N1);
            this.mSin[i] = -Math.sin(Math.PI / N1);
            N1 = N1 + N1;
        }
    }

    /**
     * Resize this {@code FFT} instance
     *
     * @param size the new length/size of this {@code FFT} instance
     */
    public void resize(int size) {
        // Make sure n is a power of 2
        if (size == 0 || (size&(size-1)) != 0)
            throw new RuntimeException("Length size must be power of 2");
        mSize = size;
        int M = (int)(Math.log(size)/Math.log(2));

        // precompute tables
        this.mCos = new double[M];
        this.mSin = new double[M];

        int N1 = 1;
        for(int i = 0; i < M; i++) {
            this.mCos[i] = Math.cos(Math.PI / N1);
            this.mSin[i] = -Math.sin(Math.PI / N1);
            N1 = N1 + N1;
        }
    }

    /**
     * Returns the Spectrum of a signal.
     *  The signal must be represented by an array.
     *  The imaginary part of the signal is assumed to be zero.
     *
     * @param  x the signal array
     *
     * @return the Frequency Spectrum
     */
    public static float[][] Spectrum(final float[] x) {
        int N = x.length;
        float[] X = new float[N];
        float[] Y = new float[N];
        System.arraycopy(x, 0, X, 0, N);
        // FFT
        fft(X, Y);
        return new float[][]{X, Y};
    }

    /**
     * Returns the Power Spectrum of a signal.
     *  The signal must be represented by an array.
     *  The imaginary part of the signal is assumed to be zero.
     *
     * @param  x the signal array
     *
     * @return the Frequency Spectrum
     */
    public static double[][] Spectrum(final double[] x) {
        int N = x.length;
        double[] X = new double[N];
        double[] Y = new double[N];
        System.arraycopy(x, 0, X, 0, N);
        // FFT
        fft(X, Y);
        return new double[][]{X, Y};
    }

    /**
     * Returns the Power Spectrum of a signal.
     *  The signal must be represented by an array.
     *  The imaginary part of the signal is assumed to be zero.
     *
     * @param  x the signal array
     *
     * @return the Power Spectrum
     */
    public static float[] PowerSpectrum(final float[] x) {
        int N = x.length;
        float[] X = new float[N];
        float[] Y = new float[N];
        System.arraycopy(x, 0, X, 0, N);
        // FFT
        fft(X, Y);
        for(int i = 0; i < N; i++){
            X[i] = (X[i]*X[i] + Y[i]*Y[i])/N;
        }
        return X;
    }

    /**
     * Returns the Power Spectrum of a signal.
     *  The signal must be represented by an array.
     *  The imaginary part of the signal is assumed to be zero.
     *
     * @param  x the signal array
     *
     * @return the Power Spectrum
     */
    public static double[] PowerSpectrum(final double[] x) {
        int N = x.length;
        double[] X = new double[N];
        double[] Y = new double[N];
        System.arraycopy(x, 0, X, 0, N);
        // FFT
        fft(X, Y);
        for(int i = 0; i < N; i++){
            X[i] = (X[i]*X[i] + Y[i]*Y[i])/N;
        }
        return X;
    }

    /**
     * Returns the Magnitude Spectrum of a signal.
     *  The signal must be represented by an array.
     *  The imaginary part of the signal is assumed to be zero.
     *
     * @param  x the signal array
     *
     * @return the Magnitude Spectrum
     */
    public static float[] MagnitudeSpectrum(final float[] x) {
        int N = x.length;
        float[] X = new float[N];
        float[] Y = new float[N];
        System.arraycopy(x, 0, X, 0, N);
        // FFT
        fft(X, Y);
        for(int i = 0; i < N; i++){
            X[i] = (float) Math.sqrt(X[i]*X[i] + Y[i]*Y[i])/N;
        }
        return X;
    }

    /**
     * Returns the Magnitude Spectrum of a signal.
     *  The signal must be represented by an array.
     *  The imaginary part of the signal is assumed to be zero.
     *
     * @param  x the signal array
     *
     * @return the Magnitude Spectrum
     */
    public static double[] MagnitudeSpectrum(final double[] x) {
        int N = x.length;
        double[] X = new double[N];
        double[] Y = new double[N];
        System.arraycopy(x, 0, X, 0, N);
        // FFT
        fft(X, Y);
        for(int i = 0; i < N; i++){
            X[i] = Math.sqrt(X[i]*X[i] + Y[i]*Y[i])/N;
        }
        return X;
    }

    /**
     * Performs the FFT on a signal.
     * The signal must be represented by a real and imaginary array of equal length.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static void fft(float[] ReX, float[] ImX){
        int N = ReX.length;
        if (N != ImX.length)
            throw new IllegalArgumentException("Real and Imaginary arrays must be identical length!");
        if ((N&(N-1)) != 0) {
            bluestein(ReX, ImX);
            return;
        }
        int M = (int)(Math.log(N)/Math.log(2));
        int ND2 = N/2;
        int J = ND2;
        float TR;
        float TI;

        // Bit Reversal Sorting
        for(int i = 1; i < N-1; i++){
            if(i < J){
                TR = ReX[J];
                TI = ImX[J];
                ReX[J] = ReX[i];
                ImX[J] = ImX[i];
                ReX[i] = TR;
                ImX[i] = TI;
            }
            int k = ND2;
            while(k <= J){
                J -= k;
                k /= 2;
            }
            J += k;
        }

        int N1; int N2 = 1;
        float SR; float SI;
        float UR; float UI;
        for(int i = 0; i < M; i++){
            N1 = N2;
            N2 = N2 + N2;
            UR = 1;
            UI = 0;
            // Calculate Sine and Cosine Values
            SR = (float) Math.cos(Math.PI/N1);
            SI = (float) -Math.sin(Math.PI/N1);
            // Loop for each Sub-DFT
            for(int j = 0; j < N1; j++){
                // Loop for each Butterfly
                for(int k = j; k < N; k+=N2){
                    int IP = k + N1;
                    // Butterfly Calculation
                    TR = ReX[IP]*UR - ImX[IP]*UI;
                    TI = ReX[IP]*UI + ImX[IP]*UR;
                    ReX[IP] = (ReX[k] - TR);
                    ImX[IP] = (ImX[k] - TI);
                    ReX[k] = (ReX[k] + TR);
                    ImX[k] = (ImX[k] + TI);
                }
                TR = UR;
                UR = TR*SR - UI*SI;
                UI = TR*SI + UI*SR;
            }
        }
    }

    /**
     * Performs the FFT on a signal.
     * The signal must be represented by a real and imaginary array of equal length.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static void fft(double[] ReX, double[] ImX){
        int N = ReX.length;
        if (N != ImX.length)
            throw new IllegalArgumentException("Real and Imaginary arrays must be identical length!");
        if ((N&(N-1)) != 0) {
            bluestein(ReX, ImX);
            return;
        }
        int M = (int)(Math.log(N)/Math.log(2));
        int ND2 = N/2;
        int J = ND2;
        double TR;
        double TI;

        // Bit Reversal Sorting
        for(int i = 1; i < N-1; i++){
            if(i < J){
                TR = ReX[J];
                TI = ImX[J];
                ReX[J] = ReX[i];
                ImX[J] = ImX[i];
                ReX[i] = TR;
                ImX[i] = TI;
            }
            int k = ND2;
            while(k <= J){
                J -= k;
                k /= 2;
            }
            J += k;
        }

        int N1; int N2 = 1;
        double SR; double SI;
        double UR; double UI;
        for(int i = 0; i < M; i++){
            N1 = N2;
            N2 = N2 + N2;
            UR = 1;
            UI = 0;
            // Calculate Sine and Cosine Values
            SR =  Math.cos(Math.PI/N1);
            SI = -Math.sin(Math.PI/N1);
            // Loop for each Sub-DFT
            for(int j = 0; j < N1; j++){
                // Loop for each Butterfly
                for(int k = j; k < N; k+=N2){
                    int IP = k + N1;
                    // Butterfly Calculation
                    TR = ReX[IP]*UR - ImX[IP]*UI;
                    TI = ReX[IP]*UI + ImX[IP]*UR;
                    ReX[IP] = (ReX[k] - TR);
                    ImX[IP] = (ImX[k] - TI);
                    ReX[k] = (ReX[k] + TR);
                    ImX[k] = (ImX[k] + TI);
                }
                TR = UR;
                UR = TR*SR - UI*SI;
                UI = TR*SI + UI*SR;
            }
        }
    }

    /**
     * Returns the FFT of the specified complex array.
     *
     * @param  x the complex array
     */
    public static void fft(Complex[] x) {
        int N = x.length;
        if ((N&(N-1)) != 0) {
            bluestein(x);
            return;
        }
        int M = (int)(Math.log(N)/Math.log(2));
        int ND2 = N/2;
        int J = ND2;
        Complex T = Complex.UNITY();

        // Bit Reversal Sorting
        for(int i = 1; i < N-1; i++){
            if(i < J){
                T = x[J];
                x[J] = x[i];
                x[i] = T;
            }
            int k = ND2;
            while(k <= J){
                J -= k;
                k /= 2;
            }
            J += k;
        }

        int N1; int N2 = 1;
        Complex S; Complex U;
        for(int i = 0; i < M; i++){
            N1 = N2;
            N2 = N2 + N2;
            U = new Complex(1, 0);
            // Calculate Sine and Cosine Values
            S = new Complex(Math.cos(Math.PI/N1), -Math.sin(Math.PI/N1));
            // Loop for each Sub-DFT
            for(int j = 0; j < N1; j++){
                // Loop for each Butterfly
                for(int k = j; k < N; k+=N2){
                    int IP = k + N1;
                    // Butterfly Calculation
                    T = x[IP].times(U);
                    x[IP] = x[k].minus(T);
                    x[k] = x[k].plus(T);
                }
                T = new Complex(U.real(), T.imag());
                U = U.times(S);
            }
        }
    }

    /**
     * Performs the Inverse FFT of a signal.
     * The signal must be represented by a real and imaginary array of equal length.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static void ifft(float[] ReX, float[] ImX) {
        fft(ImX, ReX);
        int N = ReX.length;
        for(int i = 0; i < N; i++) {
            ReX[i] /= N;
            ImX[i] /= N;
        }
    }

    /**
     * Performs the Inverse FFT of a signal.
     * The signal must be represented by a real and imaginary array of equal length.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static void ifft(double[] ReX, double[] ImX) {
        fft(ImX, ReX);
        int N = ReX.length;
        for(int i = 0; i < N; i++) {
            ReX[i] /= N;
            ImX[i] /= N;
        }
    }

    /**
     * Returns the inverse FFT of the specified complex array.
     *
     * @param  x the complex array
     * @throws IllegalArgumentException if the length of {@code x} is not a power of 2
     */
    public static void ifft(Complex[] x) {
        int n = x.length;

        // take conjugate
        for (int i = 0; i < n; i++) {
            x[i] = x[i].conj();
        }

        // compute forward FFT
        fft(x);

        // take conjugate again
        for (int i = 0; i < n; i++) {
            x[i] = x[i].conj().times(1.0/n);
        }
    }

    /**
     * Performs the FFT on a signal.
     * The signal must be represented by a 2-D array of the form {real, imaginary}.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static float[][] fftCopy(final float[] ReX, final float[] ImX){
        int N = ReX.length;
        float[] ReY = new float[N];
        float[] ImY = new float[N];
        System.arraycopy(ReX, 0 , ReY, 0, N);
        System.arraycopy(ImX, 0 , ImY, 0, N);
        fft(ReY, ImY);
        return new float[][]{ReY, ImY};
    }

    /**
     * Performs the FFT on a signal.
     * The signal must be represented by a 2-D array of the form {real, imaginary}.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static double[][] fftCopy(final double[] ReX, final double[] ImX){
        int N = ReX.length;
        double[] ReY = new double[N];
        double[] ImY = new double[N];
        System.arraycopy(ReX, 0 , ReY, 0, N);
        System.arraycopy(ImX, 0 , ImY, 0, N);
        fft(ReY, ImY);
        return new double[][]{ReY, ImY};
    }

    /**
     * Returns the FFT of the specified complex array.
     *
     * @param  x the complex array
     * @return the FFT of the complex array {@code x}
     */
    public static Complex[] fftCopy(final Complex[] x) {
        int N = x.length;
        Complex[] X = new Complex[N];
        System.arraycopy(x, 0 , X, 0, N);
        fft(X);
        return X;
    }

    /**
     * Performs the Inverse FFT on a signal.
     * The signal must be represented by a 2-D array of the form {real, imaginary}.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static float[][] ifftCopy(final float[] ReX, final float[] ImX){
        int N = ReX.length;
        float[] ReY = new float[N];
        float[] ImY = new float[N];
        System.arraycopy(ReX, 0 , ReY, 0, N);
        System.arraycopy(ImX, 0 , ImY, 0, N);
        ifft(ReY, ImY);
        return new float[][]{ReY, ImY};
    }

    /**
     * Performs the Inverse FFT on a signal.
     * The signal must be represented by a 2-D array of the form {real, imaginary}.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     * @throws IllegalArgumentException if the length of {@code ReX} is != to the length of {@code ImX}
     */
    public static double[][] ifftCopy(final double[] ReX, final double[] ImX){
        int N = ReX.length;
        double[] ReY = new double[N];
        double[] ImY = new double[N];
        System.arraycopy(ReX, 0 , ReY, 0, N);
        System.arraycopy(ImX, 0 , ImY, 0, N);
        ifft(ReY, ImY);
        return new double[][]{ReY, ImY};
    }

    /**
     * Returns the inverse FFT of the specified complex array.
     *
     * @param  x the complex array
     * @return the inverse FFT of the complex array {@code x}
     * @throws IllegalArgumentException if the length of {@code x} is not a power of 2
     */
    public static Complex[] ifftCopy(final Complex[] x) {
        int N = x.length;
        Complex[] X = new Complex[N];
        System.arraycopy(x, 0 , X, 0, N);
        ifft(X);
        return X;

    }

    /**
     * Performs the FFT on a signal using a lookup table for extra performance.
     * The signal must be represented by a real and imaginary array of equal length.
     * Requires instantiation to pre-compute the sin/cos lookup tables.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     *
     * @throws IllegalArgumentException
     *             if the length of {@code ReX} is != to the length of {@code ImX}
     *             or if the length of the arguments differs from the length of the lookup tables
     */
    public void LUfft(float[] ReX, float[] ImX) {
        if (ReX.length != ImX.length)
            throw new IllegalArgumentException("Real and Imaginary arrays must be identical length!");
        if(mSize != ReX.length)
            throw new IllegalArgumentException("This FFT object was setup to process arrays of size "
                    + mSize + " but an array of size " + ReX.length + " was given.");

        int ND2 = mSize/2;
        int M = (int)(Math.log(mSize)/Math.log(2));
        int J = ND2;
        float TR;
        float TI;

        // Bit Reversal Sorting
        for(int i = 1; i < mSize-1; i++){
            if(i < J){
                TR = ReX[J];
                TI = ImX[J];
                ReX[J] = ReX[i];
                ImX[J] = ImX[i];
                ReX[i] = TR;
                ImX[i] = TI;
            }
            int k = ND2;
            while(k <= J){
                J -= k;
                k /= 2;
            }
            J += k;
        }

        // fft
        int N1; int N2 = 1;
        float UR; float UI;
        float SR; float SI;
        for(int i = 0; i < M; i++){
            N1 = N2;
            N2 = N2 + N2;
            UR = 1;
            UI = 0;
            // Calculate Sine and Cosine Values
            SR = (float) mCos[i];
            SI = (float) mSin[i];
            // Loop for each Sub-DFT
            for(int j = 0; j < N1; j++){
                // Loop for each Butterfly
                for(int k = j; k < mSize; k+=N2){
                    int IP = k + N1;
                    // Butterfly Calculation
                    TR = ReX[IP]*UR - ImX[IP]*UI;
                    TI = ReX[IP]*UI + ImX[IP]*UR;
                    ReX[IP] = (ReX[k] - TR);
                    ImX[IP] = (ImX[k] - TI);
                    ReX[k] = (ReX[k] + TR);
                    ImX[k] = (ImX[k] + TI);
                }
                TR = UR;
                UR = TR*SR - UI*SI;
                UI = TR*SI + UI*SR;
            }
        }
    }

    /**
     * Performs the FFT on a signal using a lookup table for extra performance.
     * The signal must be represented by a real and imaginary array of equal length.
     * Requires instantiation to pre-compute the sin/cos lookup tables.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     *
     * @throws IllegalArgumentException
     *             if the length of {@code ReX} is != to the length of {@code ImX}
     *             or if the length of the arguments differs from the length of the lookup tables
     */
    public void LUfft(double[] ReX, double[] ImX) {
        if (ReX.length != ImX.length)
            throw new IllegalArgumentException("Real and Imaginary arrays must be identical length!");
        if(mSize != ReX.length)
            throw new IllegalArgumentException("This FFT object was setup to process arrays of size "
                    + mSize + " but an array of size " + ReX.length + " was given.");

        int ND2 = mSize/2;
        int M = (int)(Math.log(mSize)/Math.log(2));
        int J = ND2;
        double TR;
        double TI;

        // Bit Reversal Sorting
        for(int i = 1; i < mSize-1; i++){
            if(i < J){
                TR = ReX[J];
                TI = ImX[J];
                ReX[J] = ReX[i];
                ImX[J] = ImX[i];
                ReX[i] = TR;
                ImX[i] = TI;
            }
            int k = ND2;
            while(k <= J){
                J -= k;
                k /= 2;
            }
            J += k;
        }

        // fft
        int N1; int N2 = 1;
        double UR; double UI;
        double SR; double SI;
        for(int i = 0; i < M; i++){
            N1 = N2;
            N2 = N2 + N2;
            UR = 1;
            UI = 0;
            // Calculate Sine and Cosine Values
            SR = mCos[i];
            SI = mSin[i];
            // Loop for each Sub-DFT
            for(int j = 0; j < N1; j++){
                // Loop for each Butterfly
                for(int k = j; k < mSize; k+=N2){
                    int IP = k + N1;
                    // Butterfly Calculation
                    TR = ReX[IP]*UR - ImX[IP]*UI;
                    TI = ReX[IP]*UI + ImX[IP]*UR;
                    ReX[IP] = (ReX[k] - TR);
                    ImX[IP] = (ImX[k] - TI);
                    ReX[k] = (ReX[k] + TR);
                    ImX[k] = (ImX[k] + TI);
                }
                TR = UR;
                UR = TR*SR - UI*SI;
                UI = TR*SI + UI*SR;
            }
        }
    }

    /**
     * Returns the FFT of the specified complex array.
     *
     * @param  x the complex array
     */
    public void LUfft(Complex[] x) {
        int N = x.length;
        if ((N&(N-1)) != 0) {
            bluestein(x);
            return;
        }
        int M = (int)(Math.log(N)/Math.log(2));
        int ND2 = N/2;
        int J = ND2;
        Complex T = Complex.UNITY();

        // Bit Reversal Sorting
        for(int i = 1; i < N-1; i++){
            if(i < J){
                T = x[J];
                x[J] = x[i];
                x[i] = T;
            }
            int k = ND2;
            while(k <= J){
                J -= k;
                k /= 2;
            }
            J += k;
        }

        int N1; int N2 = 1;
        Complex S; Complex U;
        for(int i = 0; i < M; i++){
            N1 = N2;
            N2 = N2 + N2;
            U = new Complex(1, 0);
            // Calculate Sine and Cosine Values
            S = new Complex(mCos[i], mSin[i]);
            // Loop for each Sub-DFT
            for(int j = 0; j < N1; j++){
                // Loop for each Butterfly
                for(int k = j; k < N; k+=N2){
                    int IP = k + N1;
                    // Butterfly Calculation
                    T = x[IP].times(U);
                    x[IP] = x[k].minus(T);
                    x[k] = x[k].plus(T);
                }
                T = new Complex(U.real(), T.imag());
                U = U.times(S);
            }
        }
    }

    /**
     * Performs the Inverse FFT on a signal using a lookup table for extra performance.
     * The signal must be represented by a real and imaginary array of equal length.
     * Requires instantiation to pre-compute the sin/cos lookup tables.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     *
     * @throws IllegalArgumentException
     *             if the length of {@code ReX} is != to the length of {@code ImX}
     *             or if the length of the arguments differs from the length of the lookup tables
     */
    public void LUifft(float[] ReX, float[] ImX) {
        LUfft(ImX, ReX);
        int N = ReX.length;
        for(int i = 0; i < N; i++) {
            ReX[i] /= N;
            ImX[i] /= N;
        }
    }

    /**
     * Performs the Inverse FFT on a signal using a lookup table for extra performance.
     * The signal must be represented by a real and imaginary array of equal length.
     * Requires instantiation to pre-compute the sin/cos lookup tables.
     *
     * @param  ReX the real array
     * @param  ImX the imaginary array
     *
     * @throws IllegalArgumentException
     *             if the length of {@code ReX} is != to the length of {@code ImX}
     *             or if the length of the arguments differs from the length of the lookup tables
     */
    public void LUifft(double[] ReX, double[] ImX) {
        LUfft(ImX, ReX);
        int N = ReX.length;
        for(int i = 0; i < N; i++) {
            ReX[i] /= N;
            ImX[i] /= N;
        }
    }

    /**
     * Returns the inverse FFT of the specified complex array.
     *
     * @param  x the complex array
     *
     * @throws IllegalArgumentException if the length of {@code x} is not a power of 2
     */
    public void LUifft(Complex[] x) {
        int n = x.length;

        // take conjugate
        for (int i = 0; i < n; i++) {
            x[i] = x[i].conj();
        }

        // compute forward FFT
        LUfft(x);

        // take conjugate again
        for (int i = 0; i < n; i++) {
            x[i] = x[i].conj().times(1.0/n);
        }
    }

    /**
     * Computes the discrete Fourier transform (DFT) of the given complex signal, storing the result back into the vector.
     * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 fft function.
     * Uses Bluestein's chirp z-transform algorithm.
     *
     * @param  ReX real part of signal
     * @param  ImX imaginary part of signal
     * @throws IllegalArgumentException if the length of {@code ReX} does not equal
     *         the length of {@code ImX}
     */
    private static void bluestein(float[] ReX, float[] ImX) {
        // Find a power-of-2 convolution length m such that m >= n * 2 + 1
        int n = ReX.length;
        if (n != ImX.length)
            throw new IllegalArgumentException("Real and Imaginary arrays must be identical length!");
        int m = Integer.highestOneBit(n) * 4;

        // Trignometric tables
        float[] cos = new float[n];
        float[] sin = new float[n];

        // Prepare A & B convolution arrays
        float[] ReA = new float[m];
        float[] ImA = new float[m];
        float[] ReB = new float[m];
        float[] ImB = new float[m];

        ReB[0] = 1;
        ImB[0] = 0;
        for (int i = 0; i < n; i++) {
            int j = (int)((long)i * i % (n * 2));  // More accurate than j = i * i
            cos[i] = (float)Math.cos(Math.PI * j / n);
            sin[i] = (float) Math.sin(Math.PI * j / n);

            ReA[i] =  ReX[i] * cos[i] + ImX[i] * sin[i];
            ImA[i] = -ReX[i] * sin[i] + ImX[i] * cos[i];

            if(i > 0) {
                ReB[i] = ReB[m - i] = cos[i];
                ImB[i] = ImB[m - i] = sin[i];
            }
        }

        // Convolve A & B
        cconvolve(ReA, ImA, ReB, ImB);

        // Repack into X
        for (int i = 0; i < n; i++) {
            ReX[i] =  ReA[i] * cos[i] + ImA[i] * sin[i];
            ImX[i] = -ReA[i] * sin[i] + ImA[i] * cos[i];
        }
    }

    /**
	 * Computes the discrete Fourier transform (DFT) of the given complex signal, storing the result back into the vector.
	 * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 fft function.
	 * Uses Bluestein's chirp z-transform algorithm.
     *
     * @param  ReX real part of signal
     * @param  ImX imaginary part of signal
     * @throws IllegalArgumentException if the length of {@code ReX} does not equal
     *         the length of {@code ImX}
     */
    private static void bluestein(double[] ReX, double[] ImX) {
        // Find a power-of-2 convolution length m such that m >= n * 2 + 1
        int n = ReX.length;
        if (n != ImX.length)
            throw new IllegalArgumentException("Real and Imaginary arrays must be identical length!");
        int m = Integer.highestOneBit(n) * 4;

        // Trignometric tables
        double[] cos = new double[n];
        double[] sin = new double[n];

        // Prepare A & B convolution arrays
        double[] ReA = new double[m];
        double[] ImA = new double[m];
        double[] ReB = new double[m];
        double[] ImB = new double[m];

        ReB[0] = 1;
        ImB[0] = 0;
        for (int i = 0; i < n; i++) {
            int j = (int)((long)i * i % (n * 2));  // More accurate than j = i * i
            cos[i] = Math.cos(Math.PI * j / n);
            sin[i] = Math.sin(Math.PI * j/ n);

            ReA[i] =  ReX[i] * cos[i] + ImX[i] * sin[i];
            ImA[i] = -ReX[i] * sin[i] + ImX[i] * cos[i];

            if(i > 0) {
                ReB[i] = ReB[m - i] = cos[i];
                ImB[i] = ImB[m - i] = sin[i];
            }
        }

        // Convolve A & B
        cconvolve(ReA, ImA, ReB, ImB);

        // Repack into X
        for (int i = 0; i < n; i++) {
            ReX[i] =  ReA[i] * cos[i] + ImA[i] * sin[i];
            ImX[i] = -ReA[i] * sin[i] + ImA[i] * cos[i];
        }
    }

    /**
     * Computes the discrete Fourier transform (DFT) of the given complex signal, storing the result back into the vector.
     * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 fft function.
     * Uses Bluestein's chirp z-transform algorithm.
     *
     * @param  x the Complex signal array
     *
     * @throws IllegalArgumentException if the length of {@code ReX} does not equal
     *         the length of {@code ImX}
     */
    private static void bluestein(Complex[] x) {
        // Find a power-of-2 convolution length m such that M >= N * 2 + 1
        int N = x.length;
        int M = Integer.highestOneBit(N) * 4;

        // Trignometric tables
        Complex[] trig = new Complex[N];

        // Prepare A & B convolution arrays
        Complex[] A = new Complex[N];
        Complex[] B = new Complex[N];
        B[0] = Complex.UNITY();
        for (int i = 0; i < N; i++) {
            int j = (int)((long)i * i % (N * 2));  // More accurate than j = i * i
            trig[i] = new Complex(Math.cos(Math.PI * j / N), -Math.sin(Math.PI * j/ N));

            A[i] = A[i].times(trig[i]);

            if(i > 0) {
                B[i] = B[M - i] = trig[i];
            }
        }

        // Convolve A & B
        cconvolve(A, B);

        // Repack into X
        for (int i = 0; i < N; i++) {
            x[i] = A[i].times(trig[i]);
        }
    }

    /**
     * Performs the circular convolution of the two specified arrays.
     *
     * @param  ReX real part of first
     * @param  ImX imaginary part of first
     * @param  ReY real part of second
     * @param  ImY imaginary part of second
     * @throws IllegalArgumentException if the length of {@code x} does not equal
     *         the length of {@code y} or if the length is not a power of 2
     */
    public static void cconvolve(float[] ReX, float[] ImX, float[] ReY, float[] ImY) {
        // should probably pad x and y with 0s so that they have same length and are powers of 2
        int n = ReX.length;
        if (n != ReY.length)
            throw new IllegalArgumentException("Convolution arrays must be identical length!");

        // compute fft of each sequence
        fft(ReX, ImX);
        fft(ReY, ImY);

        // point-wise multiply
        for (int i = 0; i < n; i++) {
            float temp = ReX[i];
            ReX[i] = ReX[i] * ReY[i] - ImX[i] * ImY[i];
            ImX[i] = temp * ImY[i] + ImX[i] * ReY[i];
        }

        // compute inverse fft
        ifft(ReX, ImX);
    }

    /**
     * Performs the circular convolution of the two specified arrays.
     *
     * @param  ReX real part of first
     * @param  ImX imaginary part of first
     * @param  ReY real part of second
     * @param  ImY imaginary part of second
     * @throws IllegalArgumentException if the length of {@code x} does not equal
     *         the length of {@code y} or if the length is not a power of 2
     */
    public static void cconvolve(double[] ReX, double[] ImX, double[] ReY, double[] ImY) {
        // should probably pad x and y with 0s so that they have same length and are powers of 2
        int n = ReX.length;
        if (n != ReY.length)
            throw new IllegalArgumentException("Convolution arrays must be identical length!");

        // compute fft of each sequence
        fft(ReX, ImX);
        fft(ReY, ImY);

        // point-wise multiply
        for (int i = 0; i < n; i++) {
            double temp = ReX[i];
            ReX[i] = ReX[i] * ReY[i] - ImX[i] * ImY[i];
            ImX[i] = temp * ImY[i] + ImX[i] * ReY[i];
        }

        // compute inverse fft
        ifft(ReX, ImX);
    }

    /**
     * Performs the circular convolution of the two specified arrays.
     *
     * @param x the first signal array
     * @param y the second signal array
     * @throws IllegalArgumentException if the length of {@code x} does not equal
     *         the length of {@code y} or if real and imaginary lengths are not equal
     */
    public static void cconvolve(Complex[] x, Complex[] y) {
        // should probably pad x and y with 0s so that they have same length and are powers of 2
        int n = x.length;
        if (n != y.length)
            throw new IllegalArgumentException("Convolution arrays must be identical length!");

        // compute fft of each sequence
        fft(x);
        fft(y);

        // point-wise multiply
        for (int i = 0; i < n; i++) {
            x[i] = x[i].times(y[i]);
        }

        // compute inverse fft
        ifftCopy(x);
    }

    /**
     * Returns the autocorrelation of the input signal, x.
     * The imaginary part of the input signal is assumed to be zero.
     *
     * @param  x the signal array
     */
    public static float[] autocorr(final float[] x) {
        int N = x.length;
        float[] ReAC = new float[x.length*2];
        float[] ImAC = new float[x.length*2];
        System.arraycopy(x, 0, ReAC, 0, N);

        fft(ReAC, ImAC);
        for(int i = 0; i < ReAC.length; i++){
            ReAC[i] = ReAC[i]*ReAC[i] + ImAC[i]*ImAC[i];
            ImAC[i] = 0;
        }

        ifft(ReAC, ImAC);
        for(int i = 0; i < ReAC.length; i++){
            ReAC[i] /= N;
        }

        return ReAC;
    }

    /**
     * Returns the autocorrelation of the input signal, x.
     * The imaginary part of the input signal is assumed to be zero.
     *
     * @param  x the signal array
     */
    public static double[] autocorr(final double[] x) {
        int N = x.length;
        double[] ReAC = new double[x.length*2];
        double[] ImAC = new double[x.length*2];
        System.arraycopy(x, 0, ReAC, 0, N);

        fft(ReAC, ImAC);
        for(int i = 0; i < ReAC.length; i++){
            ReAC[i] = (float)(ReAC[i]*ReAC[i] + ImAC[i]*ImAC[i]);
            ImAC[i] = 0;
        }

        ifft(ReAC, ImAC);
        for(int i = 0; i < ReAC.length; i++){
            ReAC[i] /= N;
        }

        return ReAC;
    }

    /**
     * Returns the autocorrelation of the input signal, x.
     * The imaginary part of the input signal is assumed to be zero.
     *
     * @param  x the signal array
     */
    public static Complex[] autocorr(final Complex[] x) {
        int N = x.length;
        Complex[] AC = new Complex[N*2];
        float[] ReAC = new float[x.length*2];
        float[] ImAC = new float[x.length*2];
        System.arraycopy(x, 0, AC, 0, N);

        fft(AC);
        for(int i = 0; i < AC.length; i++){
            AC[i] = new Complex(AC[i].power(), 0);
        }

        ifft(AC);
        for(int i = 0; i < AC.length; i++){
            AC[i] = AC[i].divide(N);
        }

        return AC;
    }

}
