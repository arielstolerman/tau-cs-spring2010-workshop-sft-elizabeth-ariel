/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: FunctionException.java
 * Description: Code for class FunctionException, exception thrown by class Function
 * 
 */

package Function;

@SuppressWarnings("serial")
/**
 * @author Elizabeth Firman and Ariel Stolerman
 * This Exceptions are thrown by class Function on illegal construction or redefinition of N.
 */
public class FunctionException extends Exception {

	/**
	 * Constructs a FunctionException instance with the given message.
	 * @param message
	 * 			The detailed message for the exception.
	 */
	public FunctionException(String message){
		super(message);
	}

}
