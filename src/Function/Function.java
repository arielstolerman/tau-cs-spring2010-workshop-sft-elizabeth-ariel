/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: Function.java
 * Description: Code for abstract class Function, used as a parameter in class SFT
 * 
 */

package Function;

import SFT.Complex;
import SFT.FunctionException;

/**
 * This abstract class is used for describing functions over G -> C where G is Z_N1 x ... x Z_Nk or
 * alternatively G is a finite Abelian group described by (gj,Nj) j=1,...,k.
 * The SFT class uses a Function object for query access to the investigated function.
 * @author Elizabeth Firman and Ariel Stolerman
 */
public abstract class Function {
	
	// The infinity norm and Euclidean norm of the function
	protected Double infNorm;
	protected Double eucNorm;
	
	// upper bound for summation in norm calculation
	protected static double LOCAL_MAX = Double.MAX_VALUE/2.0;
	
	/* ***********
	 * constructor
	 *************/
	
	/**
	 * Constructs a Function object over G -> C for the given parameter G.
	 */
	public Function() throws FunctionException{
		this.infNorm = null;
		this.eucNorm = null;
	}
	
	/* ****************
	 * abstract methods
	 ******************/
	
	/**
	 * Returns the value of the function for the input element in G.
	 * @param elem
	 * 			The element whose this function's value is calculated.
	 * @return
	 * 			The value of the function for the input element in G.
	 */
	public abstract Complex getValue(long[] elem);
	
	/**
	 * Returns the infinity norm of this function over G.
	 * @return
	 * 			The infinity norm of this function over G.
	 */
	public abstract double calcInfinityNorm();
	
	/**
	 * Returns the Euclidean norm of this function over G.
	 * @return
	 * 			The infinity norm of this function over G.
	 */
	public abstract double calcEuclideanNorm();
}
