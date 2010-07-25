/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: SFTException.java
 * Description: Code for class SFTException, basic exception thrown by class SFT methods
 * 
 */

package SFT;

@SuppressWarnings("serial")
/**
 * This Exceptions are thrown by class SFT on illegal parameters.
 * @author Elizabeth Firman and Ariel Stolerman
 */
public class SFTException extends Exception{
	
	/**
	 * Constructs a SFTException instance with the given message.
	 * @param message
	 * 			The detailed message for the exception.
	 */
	public SFTException(String message){
		super(message);
	}
}
