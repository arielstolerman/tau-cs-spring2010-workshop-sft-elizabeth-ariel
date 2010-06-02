package Function;

import SFT.Complex;
import SFT.FunctionException;

/**
 * An abstract extension to class Function, for describing functions
 * over a finite Abelian group G &rarr; C.
 * @author Elizabeth Firman and Ariel Stolerman
 */
public abstract class FiniteAbelianFunction extends Function{
	
	// The vector of values describing an Abelian G
	protected long[][] G;
	
	/*
	 * constructors
	 */
	
	/**
	 * Constructs a Function object over G -> C for the given parameter G.
	 * @param G
	 * 			A vector of values describing G, i.e. a finite Abelian group described by (gj,Nj) j=1,...,k,
	 * 			the domain of the function. In each G[i], the first coordinate is the generator and
	 * 			the second coordinate is its order.
	 * @throws FunctionException
	 * 			If one of the given G-values is less than or equals to 0 or one of the generators is not in range.
	 */
	public FiniteAbelianFunction(long[][] G) throws FunctionException{
		for (int i=0; i<G.length; i++){
			if (G[i][1] <= 0){
				FunctionException fe = new FunctionException("all Ns must be positive.");
				throw fe;
			}
			if (G[i][0] <= -G[i][1] || G[i][0] >= G[i][1]){
				FunctionException fe = new FunctionException("all generators must be in range (-N,N).");
				throw fe;
			}
		}
		this.G = G;
		this.infNorm = null;
		this.eucNorm = null;
	}
	
	/* ****************
	 * abstract methods
	 ******************/
	
	public abstract Complex getValue(long elem);
	
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
	public long[][] getG(){
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
	public void setG(long[][] G) throws FunctionException{
		if (G.length != this.G.length)
			throw new FunctionException("the given G is of wrong length.");
		for (int i=0; i<G.length; i++){
			if (G[i][1] <= 0){
				FunctionException fe = new FunctionException("all Ns must be positive.");
				throw fe;
			}
			if (G[i][0] <= -G[i][1] || G[i][0] >= G[i][1]){
				FunctionException fe = new FunctionException("all generators must be in range (-N,N).");
				throw fe;
			}
		}
		boolean change = false;
		for (int i=0; i<G.length; i++){
			if (G[i][0] != this.G[i][0]){
				change = true;
				break;
			}
			if (G[i][1] != this.G[i][1]){
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
