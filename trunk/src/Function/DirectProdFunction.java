package Function;

import SFT.Complex;
import SFT.FunctionException;
import java.util.*;

/**
 * An abstract extension to class Function, for describing functions
 * over Z<sub>N1</sub> x ... x Z<sub>Nk</sub> &rarr; C.
 * @author Elizabeth Firman and Ariel Stolerman
 */
public abstract class DirectProdFunction extends Function{
	// The vector of values describing G, i.e. a Cartesian multiplication of Z_Ni, the domain of the function
	protected long[] G;

	/*
	 * constructors
	 */
	
	/**
	 * Constructs a function object over G -> C for the given parameter G.
	 * @param G
	 * 			A vector of values describing G, i.e. a Cartesian multiplication of Z_Ni,
	 * 			the domain of the function.
	 * @throws FunctionException
	 * 			If one of the given G-values is less than or equals to 0.
	 */
	public DirectProdFunction(long[] G) throws FunctionException{
		for (int i=0; i<G.length; i++){
			if (G[i] <= 0){
				FunctionException fe = new FunctionException("all Ns must be positive.");
				throw fe;
			}
		}
		this.G = G;
	}

	/* ****************
	 * abstract methods
	 ******************/
	
	public abstract Complex getValue(long[] elem);
	
	/* ********************
	 * non abstract methods
	 **********************/
	
	/**
	 * Returns the infinity norm of this function over G.
	 * @return
	 * 			The infinity norm of this function over G.
	 */ /*
	public double calcInfinityNorm(){
		return 1;
		//TODO
	}
	
	/**
	 * Returns the Euclidean norm of this function over G.
	 * @return
	 * 			The infinity norm of this function over G.
	 */ /*
	public double calcEuclideanNorm(){
		return 1;
		//TODO
	}
	
	/* ********************
	 * non-abstract methods
	 **********************/
	
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
		
		// run over all vectors changing the last coordinate and save the maximum norm
		for (i=0; i<currN; i++){
			currVector[coord] = i;
			if (coord == k-1){
				// recursion last step
				tmpNorm = getValue(currVector).getNorm();
			} else {
				// run recursively
				tmpNorm = calcInfinityNormRec(currVector, coord+1, k);
			}
			if (ans < tmpNorm) ans = tmpNorm;
		}
		return ans;
	}
	
	/**
	 * Returns the l2 norm of this function over G.
	 * The implementation given in the abstract class is the naive implementation.
	 * @return
	 * 			The l2 norm of this function over G.
	 */
	public double calcEuclideanNorm(){
		if (eucNorm == null){
			int k = G.length;
			double sizeOfG = 1;
			for (int i=0; i<k; i++) sizeOfG *= G[i]; 
			long[] currVector = new long[k];
			double sum = calcEuclideanNormRec(currVector, 0, k, sizeOfG);
			return Math.sqrt(sum);
		} else return eucNorm;
	}
	
	/**
	 * Recursively calculates the l2 norm of the function in a straight-forward way.
	 * @param currVector
	 * 			The vector currently checked with partial values.
	 * @param coord
	 * 			The coord to be checked for all possibilities at this iteration.
	 * @param k
	 * 			The size of G.
	 * @param sizeOfG
	 * 			The number of elements in G.
	 * @return
	 * 			The l2 norm of the function in a straight-forward way.
	 */
	private double calcEuclideanNormRec(long[] currVector, int coord, int k, double sizeOfG){
		long currN = G[coord];
		double tmpNormSquare;
		double currSum = 0;
		long i;
		
		// run over all vectors changing the last coordinate and save the sum the squre norms
		for (i=0; i<currN; i++){
			currVector[coord] = i;
			if (coord == k-1){
				// recursion last step
				tmpNormSquare = getValue(currVector).getNormSquare();
			} else {
				// run recursively
				tmpNormSquare = calcInfinityNormRec(currVector, coord+1, k);
			}
			currSum += tmpNormSquare/sizeOfG;
			if (currSum >= Double.MAX_VALUE/2) System.out.println("### attention! currSum passed max-double/2 ###");
		}
		return currSum;
	}
	
	/* *******
	 * getters
	 *********/
	
	/**
	 * Returns the vector of values describing G, the domain of the function.
	 * @return
	 * 			The vector of values describing G, the domain of the function.
	 */
	public long[] getG(){
		return G;
	}
	
	/* *******
	 * setters
	 *********/
	
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
