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
 * <p>
 *     These methods will more than likely be refactored into more specific classes in the near
 *     future.
 * </p>
 *
 * @author Marco Martinez
 */

public class Utilities {

    //TODO: Clean up and refactor

    // Do not allow instantiation
    private Utilities(){}

	/**
	 * Returns the center frequencies of the Octave-Band and Fractional-Octave-Band filter
	 * frequencies according to the
	 * <a href="https://law.resource.org/pub/us/cfr/ibr/002/ansi.s1.11.2004.pdf">
	 *     ANSI Specification
	 * </a>.
	 *
	 * @param N the fractional octave band to compute, i.e. 3 for 1/3 octave bands
	 * @param bands the number of filter bands to compute
	 *
	 * @return an array containing the center frequency of each 1/N octave band
	 */
    public static double[] ansiBands(int N, int bands) {
        double[] fi = new double[bands];
        double b = (double) 1/N;
        for(int j = 0; j < bands; j++) {
            int i = j - bands/3;
            if(N%2 == 0) fi[j] = 1000.0 * Math.pow(2.0, ((i+1)*b)/2.0);
            else fi[j] = 1000.0 * Math.pow(2.0, i*b);
        }
        return fi;
    }

	/**
	 * Returns the band limit frequencies of the Octave-Band and Fractional-Octave-Band filter
	 * frequencies according to the
	 * <a href="https://law.resource.org/pub/us/cfr/ibr/002/ansi.s1.11.2004.pdf">
	 *     ANSI Specification
	 * </a>.
	 *
	 * @param N the fractional octave band to compute, i.e. 3 for 1/3 octave bands
	 * @param bands the number of filter bands to compute
	 *
	 * @return an array of length bands+1 containing the upper and lower band limit frequency of
	 * 		   each 1/N octave band
	 */
    public static double[] ansiBandLimits(int N, int bands) {
        bands++;
        double[] fi = new double[bands];
        double b = (double) 1/N;
        for(int j = 0; j < bands; j++) {
            int i = j - bands/3;
            if(N%2 == 0) fi[j] = 710.0 * Math.pow(2.0, ((i+1)*b)/2.0);
            else fi[j] = 710.0 * Math.pow(2.0, i*b);
        }
        return fi;
    }

	/**
	 * Find the maximum of an array
	 *
	 * @param data the array
	 *
	 * @return the maximum value in the array
	 */
	public static float max(float[] data) {
        float max = -Float.MAX_VALUE;
        for (float value : data) {
            if (value > max) max = value;
        }
        return max;
    }

	/**
	 * Find the index location of the maximum value of an array
	 *
	 * @param data the array
	 *
	 * @return the index location of the maximum value in the array
	 */
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

	/**
	 * Find the average of an array
	 *
	 * @param data the array
	 *
	 * @return the average of all the values in the array
	 */
	public static float avgArray(float[] data) {
        float sum = 0;
        for (float value : data) sum += value;
        return sum/data.length;
    }

	/**
	 * Find the N highest peaks in an array
	 *
	 * @param data the array
	 * @param numPeaks the number of peaks to find
	 *
	 * @return the highest peaks in the array, in descending order
	 */
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

	/**
	 * Find the N lowest peaks in an array
	 *
	 * @param data the array
	 * @param numPeaks the number of peaks to find
	 *
	 * @return the lowest peaks in the array, in ascending order
	 */
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

	/**
	 * Find the first N peaks in an array
	 *
	 * @param data the array
	 * @param numPeaks the number of peaks to find
	 *
	 * @return the first peaks in the array, in the order they appear
	 */
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

	/**
	 * Are we currently running on Android, or regular Java?
	 *
	 * @return true if app is running on DVM (Android) and false if running on JVM (Java)
	 */
	public static boolean isAndroid(){
        return System.getProperty("java.vendor").equals("The Android Project");
    }
}
