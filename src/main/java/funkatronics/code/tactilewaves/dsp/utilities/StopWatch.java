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

/**
 * Class for timing stuff. Used in testing to measure performance and ensure real-time capability
 *
 * @author Marco Martinez
 */

public class StopWatch {

    // The time that the Stopwatch was started
    private long mStartTime = 0;
    // The time that the Stopwatch was stopped
    private long mStopTime = 0;
    // Is the Stopwatch currently running?
    private boolean mRunning = false;

    /**
     * Start the stop watch
     */
    public void start() {
        this.mStartTime = System.nanoTime();
        this.mRunning = true;
    }

    /**
     * Stop the stop watch
     */
    public void stop() {
        this.mStopTime = System.nanoTime();
        this.mRunning = false;
    }

    /**
     * Returns the amount of time (in nanoseconds) that has passed since calling {@code start()} if
     * the timer is running, or the amount of time between calls to {@code start()} and
     * {@code stop()} if the timer is not running
     *
     * @return the elapsed time (in nanoseconds)
     */
    public long getElapsedTime() {
        long elapsed;
        if (mRunning) elapsed = (System.nanoTime() - mStartTime);
        else elapsed = (mStopTime - mStartTime);
        return elapsed;
    }

    /**
     * Returns the amount of time (in seconds) that has passed since calling {@code start()} if the
     * timer is running, or the amount of time between calls to {@code start()} and {@code stop()}
     * if the timer is not running
     *
     * @return the elapsed time (in seconds)
     */
    public double getElapsedTimeSecs() {
        double elapsed;
        if (mRunning) elapsed = ((System.nanoTime() - mStartTime)/1000000000.0);
        else elapsed = ((mStopTime - mStartTime)/1000000000.0);
        return elapsed;
    }
}