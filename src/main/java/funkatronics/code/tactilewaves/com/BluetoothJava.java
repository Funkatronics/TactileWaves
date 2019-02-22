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

import javax.bluetooth.BluetoothStateException;
import javax.bluetooth.DeviceClass;
import javax.bluetooth.DiscoveryAgent;
import javax.bluetooth.DiscoveryListener;
import javax.bluetooth.LocalDevice;
import javax.bluetooth.RemoteDevice;
import javax.bluetooth.ServiceRecord;
import javax.bluetooth.UUID;
import javax.microedition.io.Connector;
import javax.microedition.io.StreamConnection;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Class for Bluetooth Serial Communication (Standard Serial COM over Bluetooth)
 * Accepting connections and managing connections are handled in separate threads to prevent any
 * blocking of the main thread.
 *
 * @author Marco Martinez
 */

public class BluetoothJava extends Bluetooth implements DiscoveryListener {

    // Declare a Logger for logging
    private static final Logger LOG = Logger.getLogger(BluetoothJava.class.getSimpleName());

    // Object used for thread synchronization
    private final static Object lock = new Object();

    // List of discovered devices
    private static Set<RemoteDevice> mDevices = new HashSet<>();

    // The URL of the RemoteDevice we wish to connect to
    private static String mURL = null;

    // The local Bluetooth device
    private LocalDevice mAdapter;

    // The Discovery Agent associated with the local device
    private DiscoveryAgent mAgent;

    // Thread to manage Bluetooth connections
    private ConnectedThread mConnectedThread;

    /**
     * Bluetooth Object Constructor
     * <p>
     *     Attempts to acquire the local Bluetooth device hardware and checks that it is correctly
     *     initialized
     * </p>
     */
    public BluetoothJava() {
        checkBluetoothState();
    }

    /**
     * Bluetooth Object Constructor
     * <p>
     *     Attempts to acquire the local Bluetooth device hardware and checks that it is correctly
     *     initialized
     * </p>
     *
     * @param listener {@link BluetoothEventListener} that is listening for Bluetooth events
     */
    public BluetoothJava(BluetoothEventListener listener) {
        mListener = listener;
        checkBluetoothState();
    }

    // Check the state of the BluetoothAdapter
    private void checkBluetoothState() {
        try {
            mAdapter = LocalDevice.getLocalDevice();
            mAgent = mAdapter.getDiscoveryAgent();
        } catch (BluetoothStateException e) {
			//LOG.warning("No Bluetooth device hardware found");
            if(mListener != null) mListener.bluetoothNotAvailable(STATE_NONE);
			setState(STATE_NONE);
            return;
        }
        if (!mAdapter.getBluetoothAddress().equals("")) {
			setState(STATE_READY);
			//LOG.info("Bluetooth device initialized: " + mAdapter.getFriendlyName() + ": " + mAdapter.getBluetoothAddress());
			System.out.println("Bluetooth device initialized: " + mAdapter.getFriendlyName() + ": " + mAdapter.getBluetoothAddress());
        }
    }

	/**
	 * Get a list of {@code RemoteDevice} objects from the computers list of paired devices
	 *
	 * @return an array of {@code RemoteDevice} objects
	 */
	private RemoteDevice[] getRemoteDevices() {
		RemoteDevice[] devices = new RemoteDevice[]{};
		if(mAgent != null)
			devices = mAgent.retrieveDevices(DiscoveryAgent.PREKNOWN);
		return devices;
	}

	public String[] getPairedDevices() {
		RemoteDevice[] devices = getRemoteDevices();
		try {
			// If there are paired devices
			if(devices != null && devices.length > 0) {
				String[] deviceList = new String[devices.length];
				// Loop through paired devices
				for (int i = 0; i < deviceList.length; i++) {
					deviceList[i] = devices[i].getFriendlyName(true) + ": "
							+ devices[i].getBluetoothAddress();
				}
				return deviceList;
			} else
				return new String[]{};
		} catch(IOException e) {
			e.printStackTrace();
		}
		return new String[]{};
	}

    public void connect(int index) {
        ConnectThread mConnectThread = new ConnectThread(getRemoteDevices()[index], this);
        mConnectThread.start();
    }

    // Manage a connected socket (StreamConnection) in a separate thread
    private void manageConnection(StreamConnection socket) {
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
        private StreamConnection mmSocket;
        private RemoteDevice mmDevice;
        private DiscoveryListener mmListener;

		/**
		 *  Construct a new {@code ConnectThread} Object
		 */
        public ConnectThread(RemoteDevice device, DiscoveryListener listener) {
            mmDevice = device;
            mmListener = listener;
        }

        public void run() {
            // Cancel discovery because it will slow down the connection
            //mAdapter.setDiscoverable(0);
            setState(STATE_CONNECTING);
            try {
                // Connect the device through the socket. This will block
                // until it succeeds or throws an exception
                mmSocket = connect();
            } catch (IOException e) {
                // Unable to connect; close the socket and get out
                e.printStackTrace();
                return;
            }

            // The connection attempt succeeded.
            // Do work to manage the connection (in a separate thread)
            manageConnection(mmSocket);
        }

		/**
		 *  Attempt to establish a SPP connection to a remote device
		 */
        private StreamConnection connect() throws IOException {
            // Get a BluetoothSocket to connect with the given BluetoothDevice
            try {
                // UUID is the app's UUID string, also used by the server code
                mAgent.searchServices(null,
                        new UUID[]{
                            new UUID("0000110100001000800000805F9B34FB", false)},
                        mmDevice,
                        mmListener);
                // Wait for service to be found....
                synchronized(lock){
                    lock.wait();
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            if(mURL == null){
                throw new IOException("Device does not support Serial Port Protocol");
            }

            return (StreamConnection)Connector.open(mURL);
        }

		/**
		 *  Cancel the attempted connection
		 */
        public void cancel() throws IOException{
            synchronized(lock){
                lock.notify();
            }
            setState(STATE_READY);
        }
    }

    /**
     * Thread to manage a connection with a bluetooth device
     */
    private class ConnectedThread extends Thread {
        private final StreamConnection mmSocket;
        private final InputStream mmInStream;
        private final OutputStream mmOutStream;

        private volatile boolean running = true;

        /**
         * Constructor
         * @param socket - socket containing connection to manage
         */
        public ConnectedThread(StreamConnection socket) {
            mmSocket = socket;
            InputStream tmpIn = null;
            OutputStream tmpOut = null;

            // Get the input and output streams, using temp objects because
            // member streams are final
            try {
                tmpIn = socket.openInputStream();
            } catch (IOException e) {
                LOG.severe("Error occurred when creating input stream: " + e.getStackTrace());
            }
            try {
                tmpOut = socket.openOutputStream();
            } catch (IOException e) {
                LOG.severe("Error occurred when creating output stream: " + e.getStackTrace());
            }

            mmInStream = tmpIn;
            mmOutStream = tmpOut;
        }

        public void run() {
            setState(STATE_CONNECTED);
            // Keep listening to the InputStream until an exception occurs
            while (running) {
                try {
                    byte[] buffer = new byte[4];
                    // Read from the InputStream
                    int byteCount = 0;
                    if(mmInStream.available() > 0) {
                        buffer = new byte[mmInStream.available()];
                        byteCount = mmInStream.read(buffer);

                        if (byteCount > 0) {
                            // Send the obtained bytes to the listener
                            sendDataToListener(buffer);
                        }
                    }
                } catch (IOException e) {
                    LOG.warning("Connection Lost" + e.getStackTrace());
                    break;
                }
            }

			try {
				cancel();
			} catch (IOException e) {
				LOG.warning("Error closing socket" + e.getStackTrace());
			}
        }

		/**
		 *  Call this to send data to the remote device
		 */
		public void write(byte aByte) throws IOException{
			if(mmOutStream != null) {
				mmOutStream.write(aByte);
			}else {
				throw new IOException("Bluetooth is not connected, no valid socket to write to.");
			}
		}

        /**
         *  Call this to shutdown the connection
         */
        public void cancel() throws IOException{
            running = false;
            mmSocket.close();
            setState(STATE_READY);
        }
    }

	//---------------------------- Interface Methods -----------------------------------------------
    public void deviceDiscovered(RemoteDevice btDevice, DeviceClass cod) {
        //add the device to the list
        if(!mDevices.contains(btDevice)){
            mDevices.add(btDevice);
        }
    }

    public void servicesDiscovered(int transID, ServiceRecord[] records) {
        if(records != null && records.length > 0){
            mURL = records[0].getConnectionURL(0,false);
        }
        synchronized (lock) {
            lock.notify();
        }
    }

    public void serviceSearchCompleted(int transID, int respCode) {
        synchronized (lock) {
            lock.notify();
        }
    }

    public void inquiryCompleted(int discType) {
        synchronized (lock) {
            lock.notify();
        }
    }

	/**
	 * Main Method used for testing and demonstration
	 *
	 * <p>
	 *     Connects to a remote device and performs handshake. Sending "hello" fro remote device
	 *     will trigger a response, while sending "quit" or "exit" from remote device will terminate
	 *     the connection. Any other messages will be printed to the screen.
	 * </p>
	 */
    public static void main(String[] args) {
        final BluetoothJava bt = new BluetoothJava(new BluetoothEventListener() {
            @Override
            public void bluetoothNotAvailable(int state) { }

            @Override
            public void bluetoothStateChanged(Bluetooth bt, int state) {
                //System.out.println("Listener state: " + state);
                if(state == Bluetooth.STATE_CONNECTED) {
                    try {
                    	System.out.println("Bluetooth connected");
                        bt.send("Connected to PC\r\n");
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                }
            }

            @Override
            public void bluetoothDataAvailable(Bluetooth bt, byte[] data) {
                String msg = new String(data);
                if(msg.toLowerCase().contains("hello")){
                    try {
                        bt.send("Hello Remote Device, from PC!\r\n");
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                    return;
                }

                if(msg.toLowerCase().contains("quit") || msg.toLowerCase().contains("exit")) {
                    try {
                        bt.terminateConnection();
                    } catch (Exception e1) {
                        e1.printStackTrace();
                    }
                    return;
                }

                try {
                    bt.send("Received: " + msg + "\n");
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }); // end of bt init

        String[] devices = bt.getPairedDevices();

		if(devices.length == 0) System.exit(0);

		System.out.println("Paired devices:");
		for(int i = 0; i < devices.length; i++)
        	System.out.println((i + 1) + ") " + devices[i]);

		Scanner keyboard = new Scanner(System.in);
		System.out.print("Select Device: ");
		int device = keyboard.nextInt() - 1;

        bt.connect(device);
    } // end of main
}
