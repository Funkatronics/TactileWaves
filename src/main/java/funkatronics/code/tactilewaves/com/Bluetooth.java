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

package funkatronics.code.tactilewaves.com;

import java.io.IOException;

/**
 * Interface that defines required Bluetooth functionality of implementing class
 * <p>
 *     Used to abstract away the different approaches for Bluetooth communication on Android and
 *     regular old Java.
 * </p>
 *
 * @author Marco Martinez
 */

public abstract class Bluetooth {

    /**
     * Bluetooth is not available or not supported
     */
    public static int STATE_NONE = -1;

    /**
     * Bluetooth is supported, but not ready/available
     */
	public static int STATE_OK = 0;

    /**
     * Bluetooth is ready to do Bluetooth stuff
     */
	public static int STATE_READY = 1;

    /**
     * Bluetooth is currently discovering
     */
	public static int STATE_DISCOVERY = 2;

    /**
     * Bluetooth is currently connecting to a device
     */
	public static int STATE_CONNECTING = 3;

    /**
     * Bluetooth is currently connected to a device
     */
	public static int STATE_CONNECTED = 4;

	// Store current state (used with state constants above)
	protected int mState;

	// Object that listens to Bluetooth events
	protected BluetoothEventListener mListener;

    /**
     * Set the {@link BluetoothEventListener} for this {@code Bluetooth} object
     *
     * @param listener a {@code BluetoothEventListener} that listens for events from this Bluetooth
     *                 object
     */
	public void setListener(BluetoothEventListener listener) {
		mListener = listener;
	}

	/**
	 * Set the Bluetooth State (INTERNAL USE ONLY)
	 *
	 * @param state STATE Constant
	 */
	protected void setState(int state) {
		if(state != mState) {
			mState = state;
			if(mListener != null) mListener.bluetoothStateChanged(this, state);
		}
	}

    /**
     * Get the current Bluetooth State
     *
     * @return integer representing the BT State
     */
	public int getState() {
		return mState;
	}

	/**
	 * Get a list of "name: address" strings from the device's list of paired devices
	 *
	 * @return an array of {@code Strings} listing the name and address of each paired device
	 */
	public abstract String[] getPairedDevices();

    /**
     * Attempt to connect to a remote bluetooth device
     *
     * @param index the index in the paired device array of the device to connect to.
     */
	public abstract void connect(int index);

    /**
     * Send a string of char's over BT (UTF-8 Encoding)
     *
     * @param send - String to send over BT
     *
     * @throws IOException if an exception is thrown when attempting to write to the stream, or if
     *                     Bluetooth is not enabled/available
     */
	public void send(String send) throws IOException {
    	send(send.getBytes("UTF-8"));
	}

    /**
	 * Send an array of bytes over BT
	 *
	 * @param send - byte array to send
	 *
	 * @throws IOException if an exception is thrown when attempting to write to the stream, or if
	 *                     Bluetooth is not enabled/available
	 */
	public void send(byte[] send) throws IOException {
		for(byte b : send)
			send(b);
	}

	/**
	 * Send a byte over BT
	 *
	 * @param send - byte to send
	 *
	 * @throws IOException if an exception is thrown when attempting to write to the stream, or if
	 *                     Bluetooth is not enabled/available
	 */
	public abstract void send(byte send) throws IOException;

    /**
     * Disconect from current Connected Thread
     *
     * @throws IOException if an exception is thrown when attempting to close the stream
     */
	public abstract void terminateConnection() throws IOException;

	// Send received data to the listener
	protected void sendDataToListener(byte[] data) {
		if(mListener != null) mListener.bluetoothDataAvailable(this, data);
	}
}
