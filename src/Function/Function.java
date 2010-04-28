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
 * @author Elizabeth Firman and Ariel Stolerman
 * This abstract class is used for describing functions over Z_N -> C. The SFT class uses a Function object
 * for query access to the investigated function.
 */
public abstract class Function {
	
	// The value describing the group Z_N, the domain of the function
	protected long N;
	// The infinity norm and Euclidean norm of the function
	private Double infNorm;
	private Double eucNorm;
	
	/*
	 * constructors
	 */
	
	/**
	 * Constructs a Function object over Z_N -> C for the given parameter N.
	 * @param N
	 * 			The value describing the group Z_N, the domain of the function
	 * @throws FunctionException
	 * 			If the given N is less than or equals to 0.
	 */
	public Function(long N) throws FunctionException{
		if (N <= 0){
			FunctionException fe = new FunctionException("N must be positive.");
			throw fe;
		}
		this.N = N;
		this.infNorm = null;
		this.eucNorm = null;
	}
	
	/*
	 * abstract methods
	 */
	
	/**
	 * Returns the value of the function for the input element in Z_N.
	 * @param elem
	 * 			The element whose this function's value is calculated.
	 * @return
	 * 			The value of the function for the input element in Z_N.
	 */
	public abstract Complex getValue(long elem);
	
	
	/*
	 * non-abstract methods
	 */
	
	/**
	 * Returns the infinity norm of this function over Z_N.
	 * The implementation given in the abstract class is the naive implementation.
	 * @return
	 * 			The infinity norm of this function over Z_N.
	 */
	public double calcInfinityNorm(){
		if (infNorm == null){
			double ans = 0;
			
			for(long i=0; i<N; i++){
				Complex val = getValue(i);
				double tmp = Math.pow(val.getRe(),2) + Math.pow(val.getIm(),2);
				if (tmp > ans) ans = tmp;
			}
			
			infNorm = Math.sqrt(ans);
		}
		return infNorm;
	}
	
	/**
	 * Returns the Euclidean norm of this function over Z_N.
	 * The implementation given in the abstract class is the naive implementation.
	 * @return
	 * 			The Euclidean norm of this function over Z_N.
	 */
	public double calcEuclideanNorm(){
		if (eucNorm == null){
			double ans = 0;
			
			for(long i=0; i<N; i++){
				Complex val = getValue(i);
				ans += (Math.pow(val.getRe(),2) + Math.pow(val.getIm(),2))/((double)N);
			}
			
			eucNorm = Math.sqrt(ans);
		}
		return eucNorm;
	}

	/*
	 * getters
	 */
	
	/**
	 * Returns the value of N, describing Z_N the domain of the function.
	 * @return
	 * 			The value of N, describing Z_N the domain of the function.
	 */
	public long getN(){
		return N;
	}
	
	/*
	 * setters
	 */
	
	/**
	 * Sets N to a new value.
	 * @param N
	 * 			The value describing the group Z_N, the domain of the function
	 * @throws FunctionException
	 * 			If the given N is less than or equals to 0.
	 */
	public void setN(long N) throws FunctionException{
		if (N <= 0){
			FunctionException fe = new FunctionException("N must be positive.");
			throw fe;
		}
		if (N != this.N){
			this.N = N;
			this.infNorm = null;
			this.eucNorm = null;
		}
	}
}
