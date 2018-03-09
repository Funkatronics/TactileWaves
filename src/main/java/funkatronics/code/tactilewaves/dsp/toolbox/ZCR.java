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

import funkatronics.code.tactilewaves.dsp.WaveFrame;

/**
 * Class to compute the rate of zero-crossings in a signal buffer
 *
 * @author Marco Martinez
 */

public class ZCR {

    /**
     * Compute the zero-crossing rateof the {@link WaveFrame}
     *
     * @param frame the {@link WaveFrame} to process
     * @return the zero-crossing rate
     */
    public static double getZCR(final WaveFrame frame){
        float[] buffer = frame.getSamples();
        int numberOfZeroCrossings = 0;
        for(int i = 1 ; i < buffer.length ; i++){
            if(buffer[i] * buffer[i-1] < 0){
                numberOfZeroCrossings++;
            }
        }

        return numberOfZeroCrossings / (double)buffer.length;
    }

    /**
     * Compute the zero-crossing rate of the signal.
     *
     * @param buffer the signal buffer to process (float array)
     * @return the zero-crossing rate
     */
    public static double getZCR(final float[] buffer){
        int numberOfZeroCrossings = 0;
        for(int i = 1 ; i < buffer.length ; i++){
            if(buffer[i] * buffer[i-1] < 0){
                numberOfZeroCrossings++;
            }
        }

        return numberOfZeroCrossings / (double)buffer.length;
    }

    /**
     * Compute the zero-crossing rate of the signal.
     *
     * @param buffer the signal buffer to process (double array)
     * @return the zero-crossing rate
     */
    public static double getZCR(final double[] buffer){
        int numberOfZeroCrossings = 0;
        for(int i = 1 ; i < buffer.length ; i++){
            if(buffer[i] * buffer[i-1] < 0){
                numberOfZeroCrossings++;
            }
        }

        return numberOfZeroCrossings / (double)buffer.length;
    }
}
