package Function;

import SFT.Complex;
import SFT.FunctionException;

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
	 */
	public double calcInfinityNorm(){
		return 1;
		//TODO
	}
	
	/**
	 * Returns the Euclidean norm of this function over G.
	 * @return
	 * 			The infinity norm of this function over G.
	 */
	public double calcEuclideanNorm(){
		return 1;
		//TODO
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
