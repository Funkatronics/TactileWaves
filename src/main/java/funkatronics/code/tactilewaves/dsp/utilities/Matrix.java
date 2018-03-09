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

import java.util.Locale;

/**
 * Class that represents Matrices. The matrix
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>12</sub></td>
 * <td>a<sub>13</sub></td>
 * <td>a<sub>14</sub></td>
 * </tr>
 * <tr>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>23</sub></td>
 * <td>a<sub>24</sub></td>
 * </tr>
 * <tr>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>34</sub></td>
 * </tr>
 * <tr>
 * <td>a<sub>41</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>44</sub></td>
 * </tr>
 * </table>
 * </p>
 * <p>
 * is stored column major in a single array, as follows:
 * </p>
 * <p>
 * <table border="1">
 * <tr>
 * <td>a<sub>11</sub></td>
 * <td>a<sub>21</sub></td>
 * <td>a<sub>31</sub></td>
 * <td>a<sub>41</sub></td>
 * <td>a<sub>12</sub></td>
 * <td>a<sub>22</sub></td>
 * <td>a<sub>32</sub></td>
 * <td>a<sub>42</sub></td>
 * <td>a<sub>13</sub></td>
 * <td>a<sub>23</sub></td>
 * <td>a<sub>33</sub></td>
 * <td>a<sub>43</sub></td>
 * <td>a<sub>14</sub></td>
 * <td>a<sub>24</sub></td>
 * <td>a<sub>34</sub></td>
 * <td>a<sub>44</sub></td>
 * </tr>
 * </table>
 * </p>
 *
 * @author Marco Martinez
 */

public class Matrix {

    // The number of rows and columns in the Matrix
    private int mRows, mCols;

    // The data within the Matrix, stored in column major order
    private double[] mData;

    /**
     * Constructs a new Matrix of specified size with all zeros
     *
     * @param rows the number of rows
     * @param cols the number of columns
     */
    public Matrix(int rows, int cols) {
        mRows = rows;
        mCols = cols;
        mData = new double[rows*cols];
    }

    /**
     * Constructs a new Matrix of specified size containing the provided data
     *
     * @param rows the number of rows
     * @param cols the number of columns
     * @param data the data to initialize the Matrix with
     * @param deep if true, a deep copy is performed, otherwise a shallow copy is made
     */
    public Matrix(int rows, int cols, double[] data, boolean deep) {
        if(rows*cols != data.length) throw new IllegalArgumentException("the length of data must equal rows*cols");
        if(deep) {
            this.mRows = rows;
            this.mCols = cols;
            this.mData = new double[rows*cols];
            System.arraycopy(data, 0, this.mData, 0, data.length);
        } else {
            this.mRows = rows;
            this.mCols = cols;
            this.mData = data;
        }
    }

    /**
     * Constructs a new Matrix of specified size containing the provided data
     *
     * @param rows the number of rows
     * @param cols the number of columns
     * @param data the data to initialize the Matrix with
     */
    public Matrix(int rows, int cols, float[] data) {
        if(rows*cols != data.length) throw new IllegalArgumentException("the length of data must equal rows*cols");
        this.mRows = rows;
        this.mCols = cols;
        this.mData = new double[rows*cols];
        for(int i = 0; i < data.length; i++)
            this.mData[i] = data[i];
    }

    /**
     * Constructs a new Matrix object from a Java 2D array
     *
     * @param data a 2D array that holds Matrix data
     */
    public Matrix(double[][] data) {
        this(data.length, data[0].length);
        for (int i = 0; i < mRows; i++)
            for (int j = 0; j < mCols; j++)
                this.mData[j* mRows + i] = data[i][j];
    }

    /**
     * Constructs a new Matrix object from a Java 2D array
     *
     * @param data a 2D array that holds Matrix data
     */
    public Matrix(float[][] data) {
        this(data.length, data[0].length);
        for (int i = 0; i < mRows; i++)
            for (int j = 0; j < mCols; j++)
                this.mData[j* mCols + i] = data[i][j];
    }

    /**
     * Constructs a new Matrix object that is a copy of another Matrix object
     *
     * @param A the Matrix to copy from
     * @param deep if true, a deep copy is performed, otherwise a shallow copy is made
     */
    public Matrix(Matrix A, boolean deep) {
        this(A.mRows, A.mCols, A.mData, deep);
    }

    /**
     * Return a N-by-N Identity Matrix
     *
     * @param N the size of the Identity Matrix
     * @return a N-by-N Identity Matrix
     */
    public static Matrix identity(int N) {
        Matrix I = new Matrix(N, N);
        for (int i = 0; i < N; i++)
            I.set(i, i, 1);
        return I;
    }

    /**
     * Get the number of rows in this matrix
     *
     * @return the number of rows
     */
    public int getRows() {
        return this.mRows;
    }

    /**
     * Get the number of columns in this matrix
     *
     * @return the number of columns
     */
    public int getCols() {
        return this.mCols;
    }

    /**
     * Get the underlying data array
     *
     * @return an array containing the Matrix values in column-major order
     */
    public double[] getData() {
        return mData;
    }

    /**
     * Get the entry at (row, col) from the Matrix
     * @param row the row location
     * @param col the column location
     * @return A(row, col)
     */
    public double get(int row, int col) {
        return mData[col* mRows + row];
    }

    /**
     * Set an entry in the Matrix to a specific value
     *
     * @param row the row location
     * @param col the column location
     * @param value the value to place, A(row, col) = value
     */
    public void set(int row, int col, double value) {
        mData[col* mRows + row] = value;
    }

    /**
     * Adds a specific value to an entry in the Matrix
     *
     * @param row the row location
     * @param col the column location
     * @param value the value to add, A(row, col) += value
     */
    public void add(int row, int col, double value) {
        mData[col* mRows + row] += value;
    }

    /**
     * Transpose this Matrix, overwriting the original
     */
    public void transpose() {
        double temp;
        for (int i = 0; i < mRows; i++)
            for (int j = 0; j < mCols; j++) {
                if(j > i) {
                    temp = this.get(i, j);
                    this.set(i, j, this.get(j, i));
                    this.set(j, i, temp);
                }
            }
    }

    /**
     * Take the inverse of this Matrix and return it
     *
     * @return the inverse of this Matrix
     */
    public Matrix inverse() {
        if (this.mCols != this.mRows) throw new RuntimeException("Matrix must be square");
        int N = this.mCols;
        Matrix I = Matrix.identity(N);
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                I.set(i, j, 0 - this.get(i, j));
            }
        }
        return I;
    }

    // swap rows a and b
    private void swap(int a, int b) {
        for(int j = 0; j < mCols; j++) {
            double temp = mData[j* mCols + a];
            mData[j* mCols + a] = mData[j* mCols + b];
            mData[j* mCols + b] = temp;
        }
    }

    /**
     * Multiply this Matrix by another Matrix
     *
     * @param B the other Matrix
     * @return the result of A*B in a new Matrix
     */
    public Matrix multiply(Matrix B) {
        Matrix A = new Matrix(this, true);
        if (A.mCols != B.mRows) throw new IllegalArgumentException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.mRows, B.mCols);
        for (int i = 0; i < A.mRows; i++)
            for (int j = 0; j < B.mCols; j++)
                for (int k = 0; k < A.mCols; k++)
                    C.add(i, j, A.get(i, k) * B.get(k, j));
        return C;
    }

    /**
     * Multiply this Matrix by a Vector
     *
     * @param x the array of values (vector)
     * @return the result of A*x in a new Matrix
     */
    public double[] multiply(double[] x) {
        Matrix A = new Matrix(this, true);
        if (A.mCols != x.length) throw new IllegalArgumentException("Illegal matrix dimensions.");
        //Matrix C = new Matrix(A.mRows, 1);
        double[] C = new double[A.mRows];
        for (int i = 0; i < A.mRows; i++)
            for (int k = 0; k < A.mCols; k++)
                C[i] += A.get(i, k) * x[k];
        return C;
    }

    /**
     * Multiply this Matrix by a Vector
     *
     * @param x the array of values (vector)
     * @return the result of A*x in a new Matrix
     */
    public Matrix multiply(float[] x) {
        Matrix A = new Matrix(this, true);
        if (A.mCols != x.length) throw new IllegalArgumentException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.mRows, 1);
        for (int i = 0; i < A.mRows; i++)
                for (int k = 0; k < A.mCols; k++)
                    C.add(i, 0, A.get(i, k) * x[k]);
        return C;
    }

    /**
     * Solve X = A\B assuming A is square and has full rank
     *
     * @param B a Matrix with the same number of rows as this Matrix
     * @return a Matrix X that has the same number of rows as this Matrix and the same number of
     *         columns as B
     */
    public Matrix solve(Matrix B) {
        if (mRows != mCols || B.mRows != mCols || B.mCols != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this, true);
        Matrix b = new Matrix(B, true);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < mCols; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < mCols; j++)
                if (Math.abs(A.get(j, i)) > Math.abs(A.get(max, i)))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular
            if (A.get(i, i) == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < mCols; j++)
                b.add(j, 0, -b.get(i, 0) * A.get(j, i) / A.get(i, i));

            // pivot within A
            for (int j = i + 1; j < mCols; j++) {
                double m = A.get(j, i) / A.get(i, i);
                for (int k = i+1; k < mCols; k++) {
                    A.add(j, k,-A.get(i, k) * m);
                }
                A.set(j, i,0.0);
            }
        }

        // back substitution
        Matrix x = new Matrix(mCols, 1);
        for (int j = mCols - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < mCols; k++)
                t += A.get(j, k) * x.get(k, 0);
            x.set(j, 0, (b.get(j, 0) - t) / A.get(j, j));
        }
        return x;

    }

    /**
     * Solve x = A\b assuming A is square and has full rank
     *
     * @param b a vector with length equal to the number of rows in this Matrix
     * @return a vector x with length equal to the number of columns in this Matrix
     */
    public double[] solve(double[] b) {
        if (mRows != mCols || b.length != mCols)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this, true);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < mCols; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < mCols; j++)
                if (Math.abs(A.get(j, i)) > Math.abs(A.get(max, i)))
                    max = j;
            A.swap(i, max);
            double temp = b[i];
            b[i] = b[max];
            b[max] = temp;

            // singular
            if (A.get(i, i) == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < mCols; j++)
                b[j] -= b[i] * A.get(j, i) / A.get(i, i);

            // pivot within A
            for (int j = i + 1; j < mCols; j++) {
                double m = A.get(j, i) / A.get(i, i);
                for (int k = i+1; k < mCols; k++) {
                    A.add(j, k,-A.get(i, k) * m);
                }
                A.set(j, i, 0.0);
            }
        }

        // back substitution
        double[] x = new double[mCols];
        for (int j = mCols - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < mCols; k++)
                t += A.get(j, k) * x[k];
            x[j] = (b[j] - t) / A.get(j, j);
        }
        return x;
    }

    /**
     * Solve x = A\b assuming A is square and has full rank
     *
     * @param b a vector with length equal to the number of rows in this Matrix
     * @return a vector x with length equal to the number of columns in this Matrix
     */
    public double[] solve(float[] b) {
        if (mRows != mCols || b.length != mCols)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this, true);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < mCols; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < mCols; j++)
                if (Math.abs(A.get(j, i)) > Math.abs(A.get(max, i)))
                    max = j;
            A.swap(i, max);
            float temp = b[i];
            b[i] = b[max];
            b[max] = temp;

            // singular
            if (A.get(i, i) == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < mCols; j++)
                b[j] -= b[i] * A.get(j, i) / A.get(i, i);

            // pivot within A
            for (int j = i + 1; j < mCols; j++) {
                double m = A.get(j, i) / A.get(i, i);
                for (int k = i+1; k < mCols; k++) {
                    A.add(j, k,-A.get(i, k) * m);
                }
                A.set(j, i, 0.0);
            }
        }

        // back substitution
        double[] x = new double[mCols];
        for (int j = mCols - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < mCols; k++)
                t += A.get(j, k) * x[k];
            x[j] = (b[j] - t) / A.get(j, j);
        }
        return x;
    }

    /**
     * Return a string representation of this Matrix
     * @return a printout of the Matrix
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        String s = "";
        for (int i = 0; i < mRows; i++) {
            for (int j = 0; j < mCols; j++)
                sb.append(String.format(Locale.getDefault(),"%2.0f ", get(i, j)));
            sb.append("\n");
        }
        return sb.toString();
    }
}
