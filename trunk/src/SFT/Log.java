package SFT;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

public class Log {
	private static String LOG_FILE = "SFT_log.txt";
	private static BufferedWriter outputFile = null;
	public static boolean IS_LOGGED = true;

	//Time
	private static final String DATE_FORMAT = "HH:mm:ss";
	private static SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT);
	private static Calendar cal = null;

	public enum DebugOutput {
		STDOUT,
		STDERR,
		FILE
	}

	/**
	 * A method for logging.
	 * The usage is: Log.log("This message goes to file..", DebugOutput.FILE);
	 * @param message: The massage explaining the current stage/problem.
	 * @param logger: The destination to write the message (STDOUT, STDERR, FILE)
	 */
	public static void log(String message, DebugOutput logger) {
		if (IS_LOGGED){
			switch (logger) {
			case STDOUT : {
				System.out.println(message);				
				break;
			}
			case STDERR : {
				System.err.println(message);    
				break;
			}
			case FILE : {
				try {
					if (outputFile == null) {
						outputFile = new BufferedWriter(new FileWriter(LOG_FILE));
					}
					cal = Calendar.getInstance();
					outputFile.write(sdf.format(cal.getTime())+" > "+message + "\r\n");
					outputFile.flush();
				} catch (IOException e) {
					System.err.println(e);
				}
				break;
			}
			}
		}
	}


	/**
	 * calls log twice with two given loggers
	 * @param message
	 * @param logger1
	 * @param logger2
	 */
	public static void log(String message, DebugOutput logger1, DebugOutput logger2) {
		log(message, logger1);
		log(message, logger2);
	}

	/**
	 * Calls log with FILE and STDOUT as DebugOutputs
	 */
	public static void log(String message){
		log(message, DebugOutput.FILE, DebugOutput.STDOUT);
	}

	/**
	 * Closing file handles.
	 */
	public static void closeLog(){
		try {
			outputFile.close();
		} catch (IOException e) {}      
	}


	// log-mode getter and setter
	
	public static boolean getLogMode(){
		return IS_LOGGED;
	}
	
	public static void setLogMode(boolean mode){
		IS_LOGGED = mode;
	}
}

