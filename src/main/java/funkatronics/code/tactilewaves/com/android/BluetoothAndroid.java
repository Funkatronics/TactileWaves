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

package funkatronics.code.tactilewaves.com.android;

import android.bluetooth.BluetoothAdapter;
import android.bluetooth.BluetoothDevice;
import android.bluetooth.BluetoothSocket;
import android.content.Context;
import android.content.DialogInterface;
import android.support.v7.app.AlertDialog;
import android.util.Log;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Set;

import funkatronics.code.tactilewaves.com.Bluetooth;
import funkatronics.code.tactilewaves.com.BluetoothEventListener;

/**
 * Class for Bluetooth Serial Communication (Standard Serial COM over Bluetooth)
 * Accepting connections and managing connections are handled in separate threads to prevent any
 * blocking of the main thread. State information is passed to the main activity for UI updates.
 *
 * @author Marco Martinez
 */

public class BluetoothAndroid extends Bluetooth {

    // Declare a TAG for logging
    private static final String TAG = BluetoothAndroid.class.getSimpleName();

    // Specific UUID for serial Bluetooth devices (RN-42, BC127, etc.)
    private java.util.UUID MY_UUID = java.util.UUID.fromString("00001101-0000-1000-8000-00805F9B34FB");

    // Bluetooth Adapter
    private BluetoothAdapter mAdapter;

	// Thread to manage a connection
    private ConnectedThread mConnectedThread;

    /**
     * Bluetooth Object Constructor
     */
    public BluetoothAndroid() {
        mAdapter = BluetoothAdapter.getDefaultAdapter();
        checkBluetoothState();
    }

    /**
     * Bluetooth Object Constructor
     *
     * @param listener {@link BluetoothEventListener} that is listening for Bluetooth events
     */
    public BluetoothAndroid(BluetoothEventListener listener) {
        mListener = listener;
        mAdapter = BluetoothAdapter.getDefaultAdapter();
        checkBluetoothState();
    }

    /**
     * Bluetooth Object Constructor
     *
     * @param listener {@link BluetoothEventListener} that is listening for Bluetooth events
     * @param adapter the {@link BluetoothAdapter} to use
     */
    public BluetoothAndroid(BluetoothEventListener listener, BluetoothAdapter adapter) {
        mListener = listener;
        //mHandler = handler;
        mAdapter = adapter;
        checkBluetoothState();
    }

    // Check the state of the BluetoothAdapter
    private void checkBluetoothState() {
        if (mAdapter == null){
            // Bluetooth Not Supported
			if(mListener != null) mListener.bluetoothNotAvailable(STATE_NONE);
            setState(STATE_NONE);
        }else{
            if (!mAdapter.isEnabled()){
				if(mListener != null) mListener.bluetoothNotAvailable(STATE_OK);
                setState(STATE_OK);
            }else{
                if(mAdapter.isDiscovering()){
                    // Wait?
                    setState(STATE_DISCOVERY);
                }else{
                    // Start BT
                    setState(STATE_READY);
                }
            }
        }
    }

	/**
	 * Get a list of {@code RemoteDevice} objects from the computers list of paired devices
	 *
	 * @return an array of {@code RemoteDevice} objects
	 */
	private BluetoothDevice[] getRemoteDevices() {
		Set<BluetoothDevice> devices = mAdapter.getBondedDevices();
		final BluetoothDevice[] devicesArray = new BluetoothDevice[devices.size()];
		devices.toArray(devicesArray);
		return devicesArray;
	}

	public String[] getPairedDevices() {
		BluetoothDevice[] devices = getRemoteDevices();
		// If there are paired devices
		if (devices != null && devices.length > 0) {
			String[] deviceList = new String[devices.length];
			// Loop through paired devices
			for (int i = 0; i < deviceList.length; i++) {
				deviceList[i] = devices[i].getName() + ": "
						+ devices[i].getAddress();
			}
			return deviceList;
		} else
			return new String[]{};
	}

    public void connect(int index) {
		ConnectThread mConnectThread = new ConnectThread(getRemoteDevices()[index]);
		mConnectThread.start();
	}

    /**
     * Choose a Bluetooth device to connect to, from the list of paired devices
     * This will open a dialog box allowing the user to choose a device to connect to out
     * of the phones paired device list
     *
	 * @param context the application context to open the dialog box within
     */
    public void chooseDeviceAndConnect(Context context) {
        if(!mAdapter.isEnabled()){
        	mListener.bluetoothNotAvailable(STATE_OK);
        }else if(getState() == STATE_READY) {

			String[] deviceList = getPairedDevices();

            AlertDialog.Builder builder = new AlertDialog.Builder(context);

            builder.setTitle("Choose Device:")
                    .setItems(deviceList, new DialogInterface.OnClickListener() {
                        public void onClick(DialogInterface dialog, int which) {
							connect(which);
                        }
                    });

            builder.setNegativeButton("Cancel", new DialogInterface.OnClickListener() {
                @Override
                public void onClick(DialogInterface dialog, int which) {
                    dialog.cancel();
                }
            }).show();
        }else {
        	mListener.bluetoothNotAvailable(mState);
        }
    }

    // Manage a connected socket in a separate thread
    private void manageConnection(BluetoothSocket socket) {
        mConnectedThread = new ConnectedThread(socket);
        mConnectedThread.start();
    }

	public void send(byte send) throws IOException {
		if(mState == STATE_CONNECTED) mConnectedThread.write(send);
	}

    public void terminateConnection() throws IOException {
        if(mState == STATE_CONNECTED){
            mConnectedThread.cancel();
            setState(STATE_READY);
        }
    }

    /**
     * Thread to create a socket to connect the bluetooth device
     * Once a successful connection is obtained, the socket is handed off
     * to a different thread (ConnectedThread) for connection management
     * NOTE: The device must first be paired
     */
    private class ConnectThread extends Thread {
        private final BluetoothSocket mmSocket;

        /**
         * Constructor
		 *
         * @param device - Bluetooth device to connect to.
         */
        public ConnectThread(BluetoothDevice device) {
            // Use a temporary object that is later assigned to mmSocket,
            // because mmSocket is final
            BluetoothSocket tmp = null;

            // Get a BluetoothSocket to connect with the given BluetoothDevice
            try {
                // MY_UUID is the app's UUID string, also used by the server code
                tmp = device.createRfcommSocketToServiceRecord(MY_UUID);
            } catch (IOException e) {
                e.printStackTrace();
            }

            mmSocket = tmp;
        }

		/**
		 * Attempt to connect to the remote device
		 */
        public void run() {
            // Cancel discovery because it will slow down the connection
            mAdapter.cancelDiscovery();
            setState(STATE_CONNECTING);
            try {
                // Connect the device through the socket. This will block
                // until it succeeds or throws an exception
                mmSocket.connect();
            } catch (IOException e) {
                // Unable to connect; close the socket and get out
                e.printStackTrace();
                cancel();
                return;
            }

            // Do work to manage the connection (in a separate thread)
            manageConnection(mmSocket);
        }

        /** Will cancel an in-progress connection, and close the socket */
        public void cancel() {
            try {
                mmSocket.close();
                setState(STATE_READY);
            } catch (IOException e) {
                Log.e(TAG, "Unable to close the client socket", e);
            }
        }
    }

    /**
     * Thread to manage a connection with a bluetooth device
     */
    private class ConnectedThread extends Thread {
        private final BluetoothSocket mmSocket;
        private final InputStream mmInStream;
        private final OutputStream mmOutStream;

        /**
         * Constructor
		 *
         * @param socket - socket containing connection to manage
         */
        public ConnectedThread(BluetoothSocket socket) {
            mmSocket = socket;
            InputStream tmpIn = null;
            OutputStream tmpOut = null;

            // Get the input and output streams, using temp objects because
            // member streams are final
            try {
                tmpIn = socket.getInputStream();
            } catch (IOException e) {
                Log.e(TAG, "Error occurred when creating input stream", e);
            }
            try {
                tmpOut = socket.getOutputStream();
            } catch (IOException e) {
                Log.e(TAG, "Error occurred when creating output stream", e);
            }

            mmInStream = tmpIn;
            mmOutStream = tmpOut;
        }

		/**
		 * Continually listen for incoming data until connection is terminated
		 */
        public void run() {
            setState(STATE_CONNECTED);
            // Keep listening to the InputStream until an exception occurs
            while (true) {
                try {
                    byte[] buffer = new byte[4];
                    // Read from the InputStream
                    int byteCount = 0;
                    if(mmInStream.available() > 0) {
                        buffer = new byte[mmInStream.available()];
                        byteCount = mmInStream.read(buffer);
                        if (byteCount > 0) {
                            Log.d(TAG, "Stream: " + new String(buffer));
                            // Send the obtained bytes to the listener
                            sendDataToListener(buffer);
                        }
                    }
					if(!isConnected()) break;
                } catch (IOException e) {
                    Log.d(TAG, "Connection Lost", e);
                    break;
                }
			}

			try {
				cancel();
			} catch (IOException e) {
				Log.w(TAG, "Error closing socket", e);
			}
		}

        /**
         * is Bluetooth currently connected?
		 *
         * @return True if connected, False otherwise
         */
        public boolean isConnected(){
            return mmSocket.isConnected();
        }

        /**
         * Call this to send data to the remote device
		 *
		 * @param aByte - byte to be sent to the remote device
         */
        public void write(byte aByte) throws IOException{
            if(mmSocket.isConnected()) {
                mmOutStream.write(aByte);
            }else {
                throw new IOException("Bluetooth is not connected, no valid socket to write to.");
            }
        }

        /**
         * Call this from the main activity to shutdown the connection
         */
        public void cancel() throws IOException{
            mmSocket.close();
            setState(STATE_READY);
        }
    }
}
