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

/**
 *  Class for windowing sample buffers
 *  <p>
 *      Arrays of both double and float are supported.
 *
 *      Hamming, Hanning, Blackman, and Gaussian Windows are available.
 *  </p>
 *
 * @author Marco Martinez
 * @see funkatronics.code.tactilewaves.dsp.WaveFrame
 */
public class Window {

    public static final int WINDOW_RECTANGULAR= -1;
    public static final int WINDOW_HAMMING = 0;
    public static final int WINDOW_HANNING = 1;
    public static final int WINDOW_BLACKMAN = 2;
    public static final int WINDOW_GAUSSIAN = 3;

    private int mWindowType = WINDOW_RECTANGULAR;

    /**
     * Construct a new {@code Window} instance with the specified window type.
     *
     * @param windowType the window type to use. Valid window types are {@code WINDOW_HAMMING},
     *                  {@code WINDOW_HANNING}, {@code WINDOW_BLACKMAN}, or {@code WINDOW_GAUSSIAN}.
     */
    public Window(int windowType) {
        this.mWindowType = windowType;
    }

    /**
     * Construct a new {@code Window} instance with a rectangular window
     */
    public Window() {
        this.mWindowType = WINDOW_RECTANGULAR;
    }

    /**
     * Get the window type assocaited wiht this {@code Window} instance.
     *
     * @return the window type. Valid window types are {@code WINDOW_HAMMING},
     *         {@code WINDOW_HANNING}, {@code WINDOW_BLACKMAN}, or {@code WINDOW_GAUSSIAN}.
     */
    public int getWindowType() {
        return mWindowType;
    }

    /**
     * Set the window type of this {@code Window} instance.
     *
     * @param windowType the window type to use. Valid window types are {@code WINDOW_HAMMING},
     *                  {@code WINDOW_HANNING}, {@code WINDOW_BLACKMAN}, or {@code WINDOW_GAUSSIAN}.
     */
    public void setWindowType(int windowType) {
        this.mWindowType = windowType;
    }

    /**
     * Sets the window function to a Hamming window for this {@code Window} instance.
     */
    public void setWindowHamming() {
        mWindowType = WINDOW_HAMMING;
    }

    /**
     * Sets the window function to a Hann window for this {@code Window} instance.
     */
    public void setWindowHanning() {
        mWindowType = WINDOW_HANNING;
    }

    /**
     * Sets the window function to a Blackman window for this {@code Window} instance.
     */
    public void setWindowBlackman() {
        mWindowType = WINDOW_BLACKMAN;
    }

    /**
     * Sets the window function to a Gaussian window for this {@code Window} instance.
     */
    public void setWindowGaussian() {
        mWindowType = WINDOW_GAUSSIAN;
    }

    /**
     * Window a signal using the window function saved to this {@code Window} instance.
     *
     * @param x the signal to e windowed
     *
     * @return the windowed signal
     */
    public float[] window(float[] x) {
        switch(mWindowType) {
            case WINDOW_HAMMING:
                return Hamming(x);
            case WINDOW_HANNING:
                return Hanning(x);
            case WINDOW_BLACKMAN:
                return Blackman(x);
            case WINDOW_GAUSSIAN:
                return Gaussian(x);
            default:
                return x;
        }
    }

    /**
     * Window a signal using the window function saved to this {@code Window} instance.
     *
     * @param x the signal to e windowed
     *
     * @return the windowed signal
     */
    public double[] window(double[] x) {
        switch(mWindowType) {
            case WINDOW_HAMMING:
                return Hamming(x);
            case WINDOW_HANNING:
                return Hanning(x);
            case WINDOW_BLACKMAN:
                return Blackman(x);
            case WINDOW_GAUSSIAN:
                return Gaussian(x);
            default:
                return x;
        }
    }

    /**
     * Returns a windowed signal using a standard Hamming window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static float[] Hamming(float[] x) {
        float[] y = new float[x.length];
        int M = x.length - 1;
        for (int i = 0; i <= M; i++) {
            y[i] = x[i];
            y[i] *= 0.54 - 0.46 * Math.cos((2 * Math.PI * i) / M);
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Hamming window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static double[] Hamming(double[] x) {
        double[] y = new double[x.length];
        int M = x.length - 1;
        for (int i = 0; i <= M; i++) {
            y[i] = x[i];
            y[i] *= 0.54 - 0.46 * Math.cos((2 * Math.PI * i) / M);
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Hann window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static float[] Hanning(float[] x) {
        float[] y = new float[x.length];
        int M = x.length - 1;
        for (int i = 0; i <= M; i++) {
            double temp = Math.sin(Math.PI * i / M);
            y[i] = x[i];
            y[i] *= (float) (temp*temp);
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Hann window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static double[] Hanning(double[] x) {
        double[] y = new double[x.length];
        int M = x.length - 1;
        for (int i = 0; i <= M; i++) {
            double temp = Math.sin(Math.PI * i / M);
            y[i] = x[i];
            y[i] *= (float) (temp*temp);
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Blackman window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static float[] Blackman(float[] x) {
        float[] y = new float[x.length];
        int M = x.length - 1;
        for (int i = 0; i <= M; i++) {
            y[i] = x[i];
            y[i] *= 0.42f - 0.5f * Math.cos((2 * Math.PI * i) / M)
                    + 0.08f * Math.cos((4 * Math.PI * i) / M);
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Blackman window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static double[] Blackman(double[] x) {
        double[] y = new double[x.length];
        int M = x.length - 1;
        for (int i = 0; i <= M; i++) {
            y[i] = x[i];
            y[i] *= 0.42f - 0.5f * Math.cos((2 * Math.PI * i) / M)
                    + 0.08f * Math.cos((4 * Math.PI * i) / M);
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Gaussian window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static float[] Gaussian(float[] x) {
        return Gaussian(x, 2.5);
    }

    /**
     * Returns a windowed signal using a standard Gaussian window.
     *
     * @param  x the signal array
     *
     * @return the windowed signal array
     */
    public static double[] Gaussian(double[] x) {
        return Gaussian(x, 2.5);
    }

    /**
     * Returns a windowed signal using a standard Gaussian window with alpha to adjust the window
     * width. Alpha is proportional to the standard deviation of the resulting window {@code sigma =
     * (N-1)/2*alpha}.
     *
     * @param x the signal array
     * @param alpha alpha value to adjust window width
     *
     * @return the windowed signal array
     */
    public static float[] Gaussian(float[] x, double alpha) {
        float[] y = new float[x.length];
        double NM1 = x.length-1;
        for (int i = 0; i <= NM1; i++) {
            double n = i - NM1/2;
            y[i] = x[i] * (float)Math.exp(-0.5*Math.pow(alpha*n/(NM1/2), 2));
        }
        return y;
    }

    /**
     * Returns a windowed signal using a standard Gaussian window with alpha to adjust the window
     * width. Alpha is proportional to the standard deviation of the resulting window {@code sigma =
     * (N-1)/2*alpha}.
     *
     * @param x the signal array
     * @param alpha alpha value to adjust window width
     *
     * @return the windowed signal array
     */
    public static double[] Gaussian(double[] x, double alpha) {
        double[] y = new double[x.length];
        double NM1 = x.length-1;
        for (int i = 0; i <= NM1; i++) {
            double n = i - NM1/2;
            y[i] = x[i] * Math.exp(-0.5*Math.pow(alpha*n/(NM1/2), 2));
        }
        return y;
    }
}
