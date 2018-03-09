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
 * Implementation of the
 * <a href="http://audition.ens.fr/adc/pdf/2002_JASA_YIN.pdf">
 *     YIN Pitch Detection Algorithm
 * </a>
 * (PDA)
 *
 *
 * @author Marco Martinez
 */

public class YIN {

    /**
     * The YIN threshold value (see paper)
     */
    private static final double mThreshold = 0.15;

    /**
     * The buffer that stores the calculated values.
     * It is exactly half the size of the input buffer.
     */
    private final float[] mYINbuffer;

    /**
     * Constructs a new YIN object to estimate the pitch (fundamental frequency) of sample buffers
     *
     * @param size the size of the sample buffers that will be passed to this YIN object
     */
    public YIN(int size) {
        mYINbuffer = new float[size/2]; //half of the buffer overlaps
    }

    /**
     * Step 2 of YIN PDA: The YINdifference function
     */
    private void YINdifference(float[] x){
        int j, tau;
        float delta;
        for(tau=0; tau < mYINbuffer.length; tau++){
            mYINbuffer[tau] = 0;
        }
        for(tau = 1; tau < mYINbuffer.length; tau++){
            for(j = 0; j < mYINbuffer.length; j++){
                delta = x[j] - x[j+tau];
                mYINbuffer[tau] += delta * delta;
            }
        }
    }

    /**
     * Step 3 of YIN PDA: The cumulative mean normalized YINdifference function
     */
    private void YINcumulativeMeanNormalizedDifference(){
        int tau;
        mYINbuffer[0] = 1;
        //start the running sum with the correct value:
        //the first value of the yinBuffer
        float runningSum = mYINbuffer[1];
        //yinBuffer[1] is always 1
        mYINbuffer[1] = 1;
        //now start at tau = 2
        for(tau = 2; tau < mYINbuffer.length; tau++){
            runningSum += mYINbuffer[tau];
            mYINbuffer[tau] *= tau / runningSum;
        }
    }

    /**
     * Step 4 of the YIN PDA
     */
    private int YINabsoluteThreshold(){
        for(int T = 1; T < mYINbuffer.length; T++){
            if(mYINbuffer[T] < mThreshold){
                while(T+1 < mYINbuffer.length && mYINbuffer[T+1] < mYINbuffer[T])
                    T++;
                return T;
            }
        }
        // no pitch found
        return -1;
    }

    /**
     * Step 5 of the YIN PDA: Refine the estimated tau value using parabolic interpolation
     * This is needed to detect higher frequencies more precisely.
     *
     * @param tauEstimate
     *            the estimated tau value.
     * @return a better, more precise tau value.
     */
    private float YINparabolicInterpolation(int tauEstimate) {
        float s0, s1, s2;
        int x0 = (tauEstimate < 1) ? tauEstimate : tauEstimate - 1;
        int x2 = (tauEstimate + 1 < mYINbuffer.length) ? tauEstimate + 1 : tauEstimate;
        if (x0 == tauEstimate)
            return (mYINbuffer[tauEstimate] <= mYINbuffer[x2]) ? tauEstimate : x2;
        if (x2 == tauEstimate)
            return (mYINbuffer[tauEstimate] <= mYINbuffer[x0]) ? tauEstimate : x0;
        s0 = mYINbuffer[x0];
        s1 = mYINbuffer[tauEstimate];
        s2 = mYINbuffer[x2];
        return tauEstimate + 0.5f * (s2 - s0 ) / (2.0f * s1 - s2 - s0);
    }

    /**
     * Step 6 of the YIN PDA
     */
    private int YINbestLocalEstimate(int tauEstimate){
        // TODO: Implement YINbestLocalEstimate as described by step 6 in the YIN paper
//        int Tmax = 25; // See YIN paper
//        float min = Float.MAX_VALUE;
//        for(int T = -Tmax/2; T <= Tmax/2; T++){
//            if(mYINbuffer[T] < min) min = mYINbuffer[T];
//        }
        tauEstimate = (int) mYINbuffer[tauEstimate];
        for(int tau = tauEstimate - tauEstimate/5; tau < tauEstimate + tauEstimate/5; tau++){
            if(mYINbuffer[tau] < mThreshold){
                while(tau+1 < mYINbuffer.length && mYINbuffer[tau+1] < mYINbuffer[tau])
                    tau++;
                return tau;
            }
        }

        return tauEstimate;
    }

    /**
     * Returns a pitch value in Hz or -1 if no pitch is detected
     *
     * @param buffer the sample buffer to analyze
     * @param sampleRate the sample rate of the audio in the buffer
     * @return a pitch value in Hz or -1 if no pitch is detected.
     */
    public float YINpitch(float[] buffer, int sampleRate) {
        int tauEstimate;
        float pitchInHertz = -1;

        //step 2
        YINdifference(buffer);

        //step 3
        YINcumulativeMeanNormalizedDifference();

        //step 4
        tauEstimate = YINabsoluteThreshold();

        //step 5
        if(tauEstimate != -1){
            float betterTau = YINparabolicInterpolation(tauEstimate);

            //step 6
            //tauEstimate = YINbestLocalEstimate(tauEstimate);

            //conversion to Hz
            pitchInHertz = sampleRate/betterTau;
        }

        return pitchInHertz;
    }
}
