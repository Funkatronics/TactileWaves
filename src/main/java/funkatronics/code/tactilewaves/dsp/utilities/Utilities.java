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

package funkatronics.code.tactilewaves.dsp.utilities;

import java.util.Arrays;

/**
 * Class that contains several utility methods such as sorting or averaging arrays
 *
 * @author Marco Martinez
 */

public class Utilities {

    //TODO: Clean up and refactor

    // Do not allow instantiation
    private Utilities(){}

    public static double[] ansiBands(int N, int chan) {
        double[] fi = new double[chan];
        double b = (double) 1/N;
        for(int j = 0; j < chan; j++) {
            int i = j - chan/3;
            if(N%2 == 0) fi[j] = 1000.0 * Math.pow(2.0, ((i+1)*b)/2.0);
            else fi[j] = 1000.0 * Math.pow(2.0, i*b);
        }
        return fi;
    }

    public static double[] ansiBandLimits(int N, int chan) {
        chan++;
        double[] fi = new double[chan];
        double b = (double) 1/N;
        for(int j = 0; j < chan; j++) {
            int i = j - chan/3;
            if(N%2 == 0) fi[j] = 710.0 * Math.pow(2.0, ((i+1)*b)/2.0);
            else fi[j] = 710.0 * Math.pow(2.0, i*b);
        }
        return fi;
    }

    public static float max(float[] aData) {
        float max = -Float.MAX_VALUE;
        for (float value : aData) {
            if (value > max) max = value;
        }
        return max;
    }

    public static int maxLoc(float[] data) {
        float max = -Float.MAX_VALUE;
        int maxLoc = 0;
        for(int i = 0; i < data.length; i++) {
            if(data[i] > max) {
                max = data[i];
                maxLoc = i;
            }
        }
        return maxLoc;
    }

    public static float avgArray(float[] data) {
        float sum = 0;
        for (float value : data) sum += value;
        return sum/data.length;
    }

    public static int[] findHighestPeaks(float[] data, int numPeaks){
        int peaks[] = new int[numPeaks];
        Arrays.fill(peaks, -1);
        float tempArray[] = new float[data.length];
        System.arraycopy(data, 0, tempArray, 0, data.length);
        int sorted[] = Sort.getSortedIndices(tempArray);
        int peaksFound = 0;
        for(int i = sorted.length - 1; i >= 0; --i){
            int maxIndex = sorted[i];
            if(maxIndex != 0 && maxIndex != sorted.length-1) {
                if ((data[maxIndex - 1] + data[maxIndex + 1]) < 2 * data[maxIndex] && data[maxIndex] > data[maxIndex-1] && data[maxIndex] > data[maxIndex+1]) {
                    peaks[peaksFound++] = maxIndex;
                    if (peaksFound >= numPeaks) break;
                }
            }
        }
        return peaks;
    }

    public static int[] findLowestPeaks(float[] data, int numPeaks){
        int peaks[] = new int[numPeaks];
        float tempArray[] = new float[data.length];
        System.arraycopy(data, 0, tempArray, 0, data.length);
        int sorted[] = Sort.getSortedIndices(tempArray);
        int peaksFound = 0;
        for (int index : sorted) {
            if (index != 0 && index != sorted.length - 1) {
                if ((data[index - 1] + data[index + 1]) < 2 * data[index] && data[index] > data[index - 1] && data[index] > data[index + 1]) {
                    peaks[peaksFound++] = index;
                    if (peaksFound >= numPeaks) break;
                }
            }
        }
        return peaks;
    }

    public static int[] findOrderedPeaks(float[] data, int numPeaks) {
        int peaks[] = new int[numPeaks];
        Arrays.fill(peaks, 0);

        int peaksFound = 0;
        for(int i = 1; i < data.length - 1; i++){
            if (data[i] > data[i-1] && data[i] > data[i+1]) {
                peaks[peaksFound++] = i;
                if (peaksFound >= numPeaks) break;
            }
        }
        return peaks;
    }

    public static boolean isAndroid(){
        return System.getProperty("java.vendor").equals("The Android Project");
    }
}
