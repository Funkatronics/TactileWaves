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
 * Class to sort arrays
 *
 * @author Marco Martinez
 */

public class Sort {

    public static int[] getSortedIndices(float[] data){
        int[] index = new int[data.length];
        for(int i = 0; i < data.length; i++)
            index[i] = i;
        sort(index, data);
        return index;
    }

    public static void sort(int[] toSort, float[] sortBy) {
        int p = 0, r = toSort.length-1;
        mergeSort(toSort, sortBy, p, r);
    }

    public static void mergeSort(int[] toSort, float[] sortBy, int p, int r) {
        if(p < r) {
            int q = (p + r) / 2;
            int[] temp = Arrays.copyOfRange(toSort, p, r + 1);
            float[] tempf = Arrays.copyOfRange(sortBy, p, r + 1);
            mergeSort(toSort, temp, sortBy, tempf, p, q);
            mergeSort(toSort, temp, sortBy, tempf, q + 1, r);
            merge(temp, toSort, tempf, sortBy, p, q, r);
        }
    }

    private static void mergeSort(int[] toSort, int[] temp, float[] sortBy, float[] tempf, int p, int r){
        if(p < r) {
            int q = (p + r) / 2;
            mergeSort(temp, toSort, tempf, sortBy, p, q);
            mergeSort(temp, toSort, tempf, sortBy,q+1, r);
            merge(toSort, temp, sortBy, tempf, p, q, r);
        }
    }

    private static void merge(int[] toSort, int[] temp, float[] sortBy, float[] tempf, int p, int q, int r) {
        int k = p;
        int L = p;
        int R = q+1;
        while(L <= q && R <= r) {
            //if(sortBy[R] > sortBy[L])
            temp[k] = sortBy[R] < sortBy[L] ? toSort[R] : toSort[L];
            tempf[k++] = sortBy[R] < sortBy[L] ? sortBy[R++] : sortBy[L++];
        }

        System.arraycopy(toSort, L, temp, k, q-L+1);
        System.arraycopy(sortBy, L, tempf, k, q-L+1);

        System.arraycopy(toSort, R, temp, k, r-R+1);
        System.arraycopy(sortBy, R, tempf, k, r-R+1);
    }
}
