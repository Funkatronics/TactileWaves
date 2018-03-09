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
 * Implements a Type-II(DCT) and Type-III(IDCT) Discrete Cosine Transform.
 * <p>
 *     Both DCT types use an FFT internally for maximized performance.
 * </p>
 *
 * @author Marco Martinez
 * @see MFCC
 */

public class DCT {

    // Do not allow instantiation
    private DCT() {}

    // Type-II DCT
    private static double[] sDCT(double[] x, boolean norm) {
        int N = x.length;
        double[] y = new double[N];
        double sum;
        for(int k = 0; k < N; k++) {
            sum = 0.0;
            for (int n = 0; n < N; n++) {
                sum += x[n] * Math.cos((2.0*(double)n + 1.0)*(double)k*Math.PI/(2.0*N));
            }
            if(norm) y[k] = sum * Math.sqrt(2.0/(double)N);
            else y[k] = sum;
        }
        if(norm) y[0] /= Math.sqrt(2.0);
        return y;
    }

    // Type-III DCT
    private static double[] sIDCT(double[] x) {
        int N = x.length;
        double[] y = new double[N];
        double sum;
        for(int k = 0; k < N; k++) {
            sum = 0.0;
            for (int n = 0; n < N; n++) {
                if(n == 0) sum += (x[n] / Math.sqrt(2.0)) * Math.cos((2.0*(double)k + 1.0)*(double)n*Math.PI/(2.0*N));
                else sum += x[n] * Math.cos((2.0*(double)k + 1.0)*(double)n*Math.PI/(2.0*N));
            }
            y[k] = sum * Math.sqrt(2.0/(double)N) ;
        }
//        y[0] /= Math.sqrt(2.0);
        return y;
    }

    /**
     * Compute a Type-II DCT of a data series
     *
     * @param x the data series array
     * @param norm should the data be normalized orthogonally?
     */
    public static void DCT(double[] x, boolean norm) {
        int N = x.length;
        double[] y = new double[N];
        for(int i = 0; i < N/2; i++){
            y[i] = x[i*2];
            y[N - i - 1] = x[i*2 + 1];
        }

        Arrays.fill(x, 0.0);
        FFT.fft(y, x);

        for(int i = 0; i < N; i++) {
            double temp = i*Math.PI/(2.0*N);
            x[i] = y[i] * Math.cos(temp) + x[i] * Math.sin(temp);
            if(norm) x[i] *= Math.sqrt(2.0/N);
        }
        if(norm) x[0] /= Math.sqrt(2.0);
    }

    /**
     * Compute a Type-III DCT (usually called the inverse DCT or IDCT) of a data series
     *
     * @param x the data series array
     * @param norm should the data be normalized orthogonally?
     */
    public static void IDCT(double[] x, boolean norm) {
        int N = x.length;
        double[] y = new double[N];
        if(N > 0) x[0] /= 2.0;
        for(int i = 0; i < N; i++) {
            if(norm) x[i] *= Math.sqrt(2.0/N);
            double temp = i*Math.PI/(2*N);
            y[i] = x[i] * Math.cos(temp);
            x[i] *= -Math.sin(temp);
            if(i == 0 && norm) y[i] *= Math.sqrt(2.0);
        }

        FFT.fft(y, x);

        for(int i = 0; i < N/2; i++){
            x[i*2] = y[i];
            x[i*2 + 1] = y[N - i - 1];
        }
        x[N - 1] = y[N/2];
    }
}
