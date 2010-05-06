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
import java.util.*;

/**
 * @author Elizabeth Firman and Ariel Stolerman
 * This abstract class is used for describing functions over G -> C where G is Z_N1 x ... x Z_Nk.
 * The SFT class uses a Function object for query access to the investigated function.
 */
public abstract class Function {
	
	// The vector of values describing G, i.e. a Cartesian multiplication of Z_Ni, the domain of the function
	protected long[] G;
	// The infinity norm and Euclidean norm of the function
	private Double infNorm;
	private Double eucNorm;
	
	// upper bound for summation in norm calculation
	private static double LOCAL_MAX = Double.MAX_VALUE/2.0;
	
	/*
	 * constructors
	 */
	
	/**
	 * Constructs a Function object over G -> C for the given parameter G.
	 * @param G
	 * 			A vector of values describing G, i.e. a Cartesian multiplication of Z_Ni,
	 * 			the domain of the function.
	 * @throws FunctionException
	 * 			If one of the given G-values is less than or equals to 0.
	 */
	public Function(long[] G) throws FunctionException{
		for (int i=0; i<G.length; i++){
			if (G[i] <= 0){
				FunctionException fe = new FunctionException("all Ns must be positive.");
				throw fe;
			}
		}
		this.G = G;
		this.infNorm = null;
		this.eucNorm = null;
	}
	
	/*
	 * abstract methods
	 */
	
	/**
	 * Returns the value of the function for the input element in G.
	 * @param elem
	 * 			The element whose this function's value is calculated.
	 * @return
	 * 			The value of the function for the input element in G.
	 */
	public abstract Complex getValue(long[] elem);
	
	
	/*
	 * non-abstract methods
	 */
	
	/**
	 * Returns the infinity norm of this function over G.
	 * The implementation given in the abstract class is the naive implementation.
	 * @return
	 * 			The infinity norm of this function over G.
	 */
	public double calcInfinityNorm(){
		if (infNorm == null){
			int k = G.length;
			long[] currVector = new long[k];
			return calcInfinityNormRec(currVector, 0, k);
		} else return infNorm;
	}
	
	/**
	 * Recursively calculates the infinity norm of the function in a straight-forward way.
	 * @param currVector
	 * 			The vector currently checked with partial values.
	 * @param coord
	 * 			The coord to be checked for all possibilities at this iteration.
	 * @param k
	 * 			The size of G.
	 * @return
	 * 			The maximum value of the function in the set of vectors with prefixes currVector.
	 */
	private double calcInfinityNormRec(long[] currVector, int coord, int k){
		long currN = G[coord];
		double ans = 0;
		double tmpNorm;
		long i;
		
		if (coord == k-1){
			// recursion last step
			// run over all vectors changing the last coordinate and save the maximum norm
			for (i=0; i<currN; i++){
				currVector[coord] = i;
				tmpNorm = getValue(currVector).getNormSquare();
				if (ans < tmpNorm) ans = tmpNorm;
			}
			ans = Math.sqrt(ans);
		} else {
			// run recursively
			for (i=0; i<currN; i++){
				currVector[coord] = i;
				tmpNorm = calcInfinityNormRec(currVector, coord+1, k);
				if (ans < tmpNorm) ans = tmpNorm;
			}
		}
		return ans;
	}
	
	
	/**
	 * Returns the Euclidean norm of this function over G.
	 * The implementation given in the abstract class is the naive implementation.
	 * @return
	 * 			The Euclidean norm of this function over G.
	 */
	public double calcEuclideanNorm(){
		if (eucNorm == null){
			int k = G.length;
			long[] currVector = new long[k];
			return Math.sqrt(calcEuclideanNormSquareRec(currVector, 0, k));
		} else return eucNorm;
	}
	
	/**
	 * Recursively calculates the Euclidean norm of the function in a straight-forward way.
	 * @param currVector
	 * 			The vector currently checked with partial values.
	 * @param coord
	 * 			The coord to be checked for all possibilities at this iteration.
	 * @param k
	 * 			The size of G.
	 * @return
	 * 			The Euclidean norm of the function in the set of vectors with prefixes currVector.
	 */
	private double calcEuclideanNormSquareRec(long[] currVector, int coord, int k){
		/*
		 * method of summation:
		 * - save partial sums in a list
		 * - for each partial sum, iterate over all Ns in G and divide
		 * - sum all results
		 * this procedure is for handling large numerical values and prevent overflow.
		 */
		
		long currN = G[coord];
		double ans = 0;
		long i;
		
		if (coord == k-1){
			// recursion last step
			// run over all vectors changing the last coordinate and calculate the next element in the sum
			List<Double> partialSums = new ArrayList<Double>();
			for (i=0; i<currN; i++){
				currVector[coord] = i;
				ans += getValue(currVector).getNormSquare();
				if (ans > LOCAL_MAX){
					partialSums.add(ans);
					ans = 0;
				}
			}
			// calculate the contribution of all partial sums for the prefix (x_1,...,x_(k-1),?) to the norm
			for(int j=0; j<partialSums.size(); i++){
				double partialSum = partialSums.remove(j);
				for(long N:G){
					partialSum /= (double)N;
				}
				partialSums.add(partialSum);
			}
			ans = 0;
			for (double partialSum: partialSums){
				ans += partialSum;
			}
		} else {
			// run recursively
			for (i=0; i<currN; i++){
				currVector[coord] = i;
				ans += calcEuclideanNormSquareRec(currVector, coord+1, k);
			}
		}
		return ans;
	}

	/*
	 * getters
	 */
	
	/**
	 * Returns the vector of values describing G, the domain of the function.
	 * @return
	 * 			The vector of values describing G, the domain of the function.
	 */
	public long[] getG(){
		return G;
	}
	
	/*
	 * setters
	 */
	
	/**
	 * Sets G to a new value.
	 * @param G
	 * 			The vector of values describing G, the domain of the function
	 * @throws FunctionException
	 * 			If the one of the given values is less than or equals to 0.
	 */
	public void setG(long[] G) throws FunctionException{
		if (G.length != this.G.length)
			throw new FunctionException("the given G is of wrong length.");
		for (int i=0; i<G.length; i++){
			if (G[i] <= 0){
				FunctionException fe = new FunctionException("all Ns must be positive.");
				throw fe;
			}
		}
		boolean change = false;
		for (int i=0; i<G.length; i++){
			if (G[i] != this.G[i]){
				change = true;
				break;
			}
		}
		if (change){
			this.G = G;
			this.infNorm = null;
			this.eucNorm = null;
		}
	}
}
