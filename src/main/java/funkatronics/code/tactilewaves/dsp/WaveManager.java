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

package funkatronics.code.tactilewaves.dsp;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

import funkatronics.code.tactilewaves.dsp.toolbox.Window;
import funkatronics.code.tactilewaves.io.WaveFormat;
import funkatronics.code.tactilewaves.io.WaveInputStream;

/**
 * Class that manages a {@link WaveInputStream} to acquire (optionally) overlapping and windowed
 * frames of audio and pass them to a chain of processors.
 * <p>
 *     A {@link WaveFrameListener} is used to send completed (processed) frames to anything that
 *     needs to listen for these events. Useful for sending processed data back to the main thread.
 * </p>
 *
 * @author Marco Martinez
 *
 * @see WaveFrame
 * @see WaveFrameListener
 */

public class WaveManager implements Runnable{

    private static final Logger LOG = Logger.getLogger(WaveManager.class.getSimpleName());

    // Is the WaveManager running?
    private boolean mRunning = false;

    // Input stream to read audio from
    private WaveInputStream mInput;

    // Window object to window audio frames
    private Window mWindow;

    // The format of the audio being processed
    private WaveFormat mFormat;

    // List of listener objects to send messages/events to
    private List<WaveFrameListener> mListeners;

    // The actual audio samples
    private float[] mSamples;

    // The amount of overlap between buffers, in samples
    private int mOverlap;
    // The length of the frame (sample buffer)
    private int mLength;

    // The number of frames that have been successfully processed by this object
    private int mFramesProcessed = 0;
    // The total number of samples that this object has read from the input stream
    private int mTotalSamplesRead = 0;

    // Effects chain - a linked list of WaveProcessor objects that execute sequentially.
    private LinkedList<WaveProcessor> mFXChain;

    /**
     * Construct a new {@code WaveManager} object from a {@link WaveInputStream} with the specified
     * frame length and overlap.
     *
     * @param input the {@link WaveInputStream} to read audio from
     * @param frameLength the length of each frame, in samples
     * @param overlap the length of the overlap of each frame, in samples
     */
    public WaveManager(WaveInputStream input, int frameLength, int overlap) {
        this.mInput = input;
        this.mFormat = input.getFormat();
        this.mLength = frameLength;
        this.mOverlap = overlap;
        this.mFramesProcessed = 0;
        this.mTotalSamplesRead = 0;
        this.mRunning = false;
        this.mWindow = new Window(Window.WINDOW_RECTANGULAR);
        this.mSamples = new float[frameLength];
        this.mFXChain = new LinkedList<>();
    }

    /**
     * Add a {@link WaveFrameListener} to the list of listeners
     * <p>
     *     A {@link WaveFrameListener} is used to send completed (processed) frames to anything that
     *     needs to listen for these events. Useful for sending processed data back to the main thread.
     * </p>
     *
     * @param listener the {@code WaveFrameListener} that should listen to events from this object
     */
    public void addListener(WaveFrameListener listener) {
        if(!mListeners.contains(listener))
            mListeners.add(listener);
    }

    /**
     * Remove a {@link WaveFrameListener} from the list of listeners
     *
     * @param listener the {@code WaveFrameListener} to remove (if it exists)
     */
    public void removeListener(WaveFrameListener listener) {
        mListeners.remove(listener);
    }

    /**
     * Add a {@link WaveProcessor} object to the processing chain
     *
     * @param effect the {@code WaveProcessor} to add to the end of the chain
     */
    public void addEffectToChain(WaveProcessor effect) {
        mFXChain.add(effect);
    }

    /**
     * Remove a {@link WaveProcessor} object from the processing chain (if it exists in the chain)
     *
     * @param effect the {@code WaveProcessor} to remove from the chain
     */
    public void removeEffectFromChain(WaveProcessor effect) {
        mFXChain.remove(effect);
        effect.processingFinished();
    }

    /**
     * Is this {@code WaveManager} currently running?
     *
     * @return true if the {@code WaveManager} is running, false otherwise
     */
    public boolean isRunning() {
        return mRunning;
    }

    /**
     * Start processing audio with this {@code WaveManager} object
     */
    public void start() {
        new Thread(this).start();
    }

    /**
     * Stop processing audio with this {@code WaveManager} object
     */
    public void stop() {
        mInput.stop();
        mRunning = false;
    }

    /**
     * Where the actual reading and processing of audio takes place.
     * <p>
     *     Calling {@code run()} directly will begin processing, but this method will block until
     *     {@code stop()} is called, so it is recommended to place every {@code WaveManager} in a
     *     dedicated thread, separate from the main thread.
     * </p>
     */
    public void run() {
        LOG.info("WaveManager starting...");
        mRunning = true;
        mInput.start();

        int samplesRead = 0;
        do {
            WaveFrame frame = new WaveFrame(mFormat);
            try {
                samplesRead = readNextFrame();
                frame.updateFrame(mWindow.window(mSamples));
            } catch (IOException e) {
                throw new Error(e);
            }
            // Send frame to each processing object
            for (WaveProcessor effect : mFXChain)
                effect.process(frame);
            // Send processed frame to all listeners
            for(WaveFrameListener listener : mListeners)
                listener.newFrameAvailable(frame);
            // Increment the number of frames processed
            mFramesProcessed++;
        } while(samplesRead != 0 && mRunning);

        LOG.info("WaveManager stopping...");
    }

    // Read the next frame of audio data from the WaveInputStream
    private int readNextFrame() throws IOException {
        // If this is Not the first frame, shift overlapping samples
        if(mTotalSamplesRead != 0) {
            System.arraycopy(mSamples, mOverlap, mSamples, 0, mOverlap);
        }

        int samplesRead = 0;
        while(samplesRead < mLength - mOverlap) {
            int read = mInput.read(mSamples, mOverlap +samplesRead, mLength -samplesRead);
            if(read == -1) {
                // Zero pad the sample buffer if the end of the stream is reached
                for(int i = mOverlap + samplesRead; i < mLength; i++)
                    mSamples[i] = 0;
                mTotalSamplesRead += samplesRead;
                return samplesRead;
            }
            else samplesRead += read;
        }

        if(samplesRead != mLength - mOverlap)
            throw new IOException("An IO error has occurred on the input stream: INVALID_STATE");

        mTotalSamplesRead += samplesRead;
        return samplesRead;
    }

}
