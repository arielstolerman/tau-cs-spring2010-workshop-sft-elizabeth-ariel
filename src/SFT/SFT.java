/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: SFT.java
 * Description: Code for class SFT, the SFT algorithm implementation
 * 
 */

package SFT;

import java.io.*;
import java.util.*;

import Function.*;
import SFT.SFTUtils.*;

/**
 * An implementation of the SFT algorithm for finding the list of elements whose Fourier
 * coefficients are significant for a given function &fnof;: G &rarr; C, where G is a Cartesian product of
 * finite groups (i.e. Z<sub>N1</sub> x ... x Z<sub>Nk</sub>) described by a list of N<sub>j</sub>'s or an
 * Abelian group described by a list of N<sub>j</sub>'s and the corresponding generators g<sub>j</sub>'s.<br>
 * <br>
 * The implementation is based on the SFT algorithm described in "Learning Noisy Characters, Multiplication Codes, 
 * And Cryptographic Hardcore Predicates" (Adi Akavia, 2008, page 52).<br>
 * <br>
 * This library is a project in CS workshop, TAU, Spring 2010.
 * @author Elizabeth Firman and Ariel Stolerman
 */
public class SFT {
	
	/* **************
	 * default values
	 * **************/
	protected static final int DEFAULT_NUM_OF_ITERATIONS = 1;
	
	/* *******************************************************************
	 * Interface public functions for direct product G (Z_N1 x ... x Z_Nk)
	 *********************************************************************/
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * confidence set by the selection of the m_A and m_B values.
	 * Here m_A and m_B are given directly by the user for calculating the sizes of groups A, Btl where
	 * Q = {x-y | x in A, y in Btl}, instead of the original algorithm that calculates m_A and m_B by given
	 * parameters.
	 * This version of the function runs only one iteration of the original SFT algorithm.
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose elements and their significant coefficients are returned.
	 * 				Used for query access.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 * @throws FunctionException
	 * 				If the creation of the difference function between iterations of the SFT procedure is invalid.
	 * 				Should not be thrown in this version.
	 */
	public static Map<long[],Complex> getSignificantElements(long[] G, double tau, DirectProdFunction func,
			long m_A, long m_B) throws SFTException,FunctionException{
		return getSignificantElements(G,tau,func,m_A,m_B,DEFAULT_NUM_OF_ITERATIONS);
	}
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * confidence set by the selection of the m_A and m_B values.
	 * Here m_A and m_B are given directly by the user for calculating the sizes of groups A, Btl where
	 * Q = {x-y | x in A, y in Btl}, instead of the original algorithm that calculates m_A and m_B by given
	 * parameters. 
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose elements and their significant coefficients are returned.
	 * 				Used for query access.
	 * @param numOfIterations
	 * 				The number of SFT procedure iterations to run. Each iteration is ran with the difference function
	 * 				of the given function and the output of the previous SFT iteration.
	 * 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
	 * 				with greater precision.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 * @throws FunctionException
	 * 				If the creation of the difference function between iterations of the SFT procedure is invalid.
	 */
	public static Map<long[],Complex> getSignificantElements(long[] G, double tau, DirectProdFunction func,
			long m_A, long m_B, int numOfIterations) throws SFTException,FunctionException{
		// skip the calculation of m_A and m_B and call the main algorithm directly
		return runMainSFTAlgorithm(G, tau, func, numOfIterations, m_A, m_B);
	}
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * delta-confidence.
	 * The algorithm includes a calculation of the error-bound, based on the delta-input and a some constant.
	 * This implementation allows the user (who knows the algorithm) to state this constant.
	 * The algorithm also includes a calculation of randomly generated sets of elements in G,
	 * of sizes defined as m_A and m_B in the paper, that uses some constant. This implementation allows the user
	 * (who knows the algorithm) to state this constant as well.
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
	 * @param numOfIterations
	 * 				The number of SFT procedure iterations to run. Each iteration is ran with the difference function
	 * 				of the given function and the output of the previous SFT iteration.
	 * 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
	 * 				with greater precision.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param fInfNorm
	 * 				The infinity norm of the function.
	 * @param fEuclideanNorm
	 * 				The Euclidean norm of the function.
	 * @param deltaCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @param randSetsCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 * @throws FunctionException
	 * 				If the creation of the difference function between iterations of the SFT procedure is invalid.
	 */
	public static Map<long[],Complex> getSignificantElements(long[] G, double tau, DirectProdFunction func,
			int numOfIterations, double delta, double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff)
			throws SFTException,FunctionException{
		// first calculate m_A and m_B as described in algorithms 3.9 and 3.10
		long[] randSetsSizes = getRandomSetsSizes(G,delta,tau,fInfNorm,fEuclideanNorm,deltaCoeff,fEuclideanNorm);
		// call main algorithm
		return runMainSFTAlgorithm(G, tau, func, numOfIterations, randSetsSizes[0], randSetsSizes[1]);
	}
	
	/* *****************************************************
	 * Interface public functions for finite Abelian G
	 *******************************************************/
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * confidence set by the selection of the m_A and m_B values.
	 * Here m_A and m_B are given directly by the user for calculating the sizes of groups A, Btl where
	 * Q = {x-y | x in A, y in Btl}, instead of the original algorithm that calculates m_A and m_B by given
	 * parameters.
	 * This version of the function runs only one iteration of the original SFT algorithm.
	 * @param G
	 * 				The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
	 * 				corresponding generators for Nj.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose elements and their significant coefficients are returned.
	 * 				Used for query access.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 * @throws FunctionException
	 * 				If the creation of the difference function between iterations of the SFT procedure is invalid.
	 * 				Should not be thrown in this version.
	 */
	public static Map<Long,Complex> getSignificantElements(long[][] G, double tau, FiniteAbelianFunction func,
			long m_A, long m_B) throws SFTException{
		return getSignificantElements(G,tau,func,DEFAULT_NUM_OF_ITERATIONS,m_A,m_B);
	}
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * confidence set by the selection of the m_A and m_B values.
	 * Here m_A and m_B are given directly by the user for calculating the sizes of groups A, Btl where
	 * Q = {x-y | x in A, y in Btl}, instead of the original algorithm that calculates m_A and m_B by given
	 * parameters. 
	 * @param G
	 * 				The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
	 * 				corresponding generators for Nj.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose elements and their significant coefficients are returned.
	 * 				Used for query access.
	 * @param numOfIterations
	 * 				The number of SFT procedure iterations to run. Each iteration is ran with the difference function
	 * 				of the given function and the output of the previous SFT iteration.
	 * 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
	 * 				with greater precision.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 * @throws FunctionException
	 * 				If the creation of the difference function between iterations of the SFT procedure is invalid.
	 */
	public static Map<Long,Complex> getSignificantElements(long[][] G, double tau, FiniteAbelianFunction func,
			int numOfIterations, long m_A, long m_B) throws SFTException{
		// create parameters for the direct product version
		long[] dpG = SFTUtils.getGFromAbelianFunc(func);
		try{
			DirectProdFunction dpFunc = new DirectedProdFromAbelianFunc(dpG,func);
			// call the direct product version of this method
			Map<long[],Complex> Ltag = getSignificantElements(dpG,tau,dpFunc,m_A,m_B,numOfIterations);
			// return the finite Abelian representation of the result set
			return SFTUtils.calcElemCoeffPairsAbelian(Ltag, G);
		} catch (FunctionException fe){
			System.err.println("Invalid function.");
			return null;
		}
	}
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * delta-confidence.
	 * The algorithm includes a calculation of the error-bound, based on the delta-input and a some constant.
	 * This implementation allows the user (who knows the algorithm) to state this constant.
	 * The algorithm also includes a calculation of randomly generated sets of elements in G,
	 * of sizes defined as m_A and m_B in the paper, that uses some constant. This implementation allows the user
	 * (who knows the algorithm) to state this constant as well.
	 * @param G
	 * 				The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
	 * 				corresponding generators for Nj.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
	 * @param numOfIterations
	 * 				The number of SFT procedure iterations to run. Each iteration is ran with the difference function
	 * 				of the given function and the output of the previous SFT iteration.
	 * 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
	 * 				with greater precision.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param fInfNorm
	 * 				The infinity norm of the function.
	 * @param fEuclideanNorm
	 * 				The Euclidean norm of the function.
	 * @param deltaCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @param randSetsCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 * @throws FunctionException
	 * 				If the creation of the difference function between iterations of the SFT procedure is invalid.
	 */
	public static Map<Long,Complex> getSignificantElements(long[][] G, double tau, FiniteAbelianFunction func,
			int numOfIterations, double delta, double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff)
			throws SFTException{
		// create parameters for the direct product version
		long[] dpG = SFTUtils.getGFromAbelianFunc(func);
		try{
			DirectProdFunction dpFunc = new DirectedProdFromAbelianFunc(dpG,func);
			// call the direct product version of this method
			Map<long[],Complex> Ltag =
				getSignificantElements(dpG,tau,dpFunc,numOfIterations,delta,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);
			// return the finite Abelian representation of the result set
			return SFTUtils.calcElemCoeffPairsAbelian(Ltag, G);
		} catch (FunctionException fe){
			System.err.println("Invalid function.");
			return null;
		}
	}	
	
	/* ************************************************************************
	 * Matlab interface public functions for Cartesian product of finite groups
	 * ************************************************************************/
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static MatlabTemporaryRepositoryDirectProd runMatlabSFTPart1Internal(Boolean isLogged, Long[] G, double tau,
			long m_A, long m_B) throws SFTException{
		// set variables to fit algorithm
		Log.setLogMode(isLogged);
		long[] g = new long[G.length];
		for(int i=0; i<G.length; i++) g[i] = G[i];
		
		// run algorithm
		Set<long[]>[][] sets = runMatlabSFTPart1Internal(g, tau, m_A, m_B);
		
		// fit results to Matlab
		Set<long[]> longQ = sets[sets.length-1][0];
		Long[][] Q = new Long[longQ.size()][];
		int i=0;
		for (long[] elem: longQ){
			int len = elem.length;
			Long[] Elem = new Long[len];
			for (int k=0; k<len; k++) Elem[k] = new Long(elem[k]);
			Q[i++] = Elem;
		}
		
		// return sets to be used by part 2, Q to be used in the matlab query calculation and null for the query
		// (will be set in the matlab script to be passed to part 2)
		MatlabTemporaryRepositoryDirectProd matlabRep = new MatlabTemporaryRepositoryDirectProd(sets, Q, null);
		return matlabRep;
	}
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static MatlabTemporaryRepositoryDirectProd runMatlabSFTPart1Internal(Boolean isLogged, Long[] G, double tau,
			double delta_t, double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff)
	throws SFTException{
		// set variables to fit algorithm
		Log.setLogMode(isLogged);
		long[] g = new long[G.length];
		for(int i=0; i<G.length; i++) g[i] = G[i];
		
		long[] randSetsSizes = getRandomSetsSizes(g,delta_t,tau,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);
		
		// run algorithm
		Set<long[]>[][] sets = runMatlabSFTPart1Internal(g, tau, randSetsSizes[0], randSetsSizes[1]);
		
		// fit results to Matlab
		Set<long[]> longQ = sets[sets.length-1][0];
		Long[][] Q = new Long[longQ.size()][];
		int i=0;
		for (long[] elem: longQ){
			int len = elem.length;
			Long[] Elem = new Long[len];
			for (int k=0; k<len; k++) Elem[k] = new Long(elem[k]);
			Q[i++] = Elem;
		}
		
		// return sets to be used by part 2, Q to be used in the matlab query calculation and null for the query
		// (will be set in the matlab script to be passed to part 2)
		MatlabTemporaryRepositoryDirectProd matlabRep = new MatlabTemporaryRepositoryDirectProd(sets, Q, null);
		return matlabRep;
	}
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static SFTUtils.MatlabTemporaryResultDirectProd runMatlabSFTPart2Internal(Long[] G, double tau,
			MatlabTemporaryRepositoryDirectProd matlabRep, int numOfIterations) throws SFTException,FunctionException{
		// fit parameters to java
		long[] g = new long[G.length];
		for(int i=0; i<G.length; i++) g[i] = G[i];
		Set<long[]>[][] sets = matlabRep.getSets();
		Map<String,Complex> query = (HashMap<String,Complex>)matlabRep.getQuery();

		// run part 2 and return the results fitted to matlab (Long[] elements)
		Map<long[],Complex> tmpRes = runMatlabSFTPart2Internal(g, tau, sets, query, numOfIterations);
		int setSize = tmpRes.keySet().size();
		Long[][] keys = new Long[setSize][];
		Complex[] values = new Complex[setSize];
		int size = G.length;
		int index = 0;
		for (long[] elem: tmpRes.keySet()){
			Long[] e = new Long[size];
			for (int i=0; i<size; i++) e[i] = new Long(elem[i]);
			
			// update results
			keys[index] = e;
			values[index] = tmpRes.get(elem);
			index++;
		}
		System.out.println("<<<<<<<<<<<<<<<<<<<< "+values[2]+" >>>>>>>>>>>>>>>>>>>>>>>");
		return new MatlabTemporaryResultDirectProd(keys, values);
	}
	
	/* ***********************************************************
	 * Matlab interface public functions for finite Abelian groups
	 * ***********************************************************/
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static MatlabTemporaryRepositoryFiniteAbelian runMatlabSFTPart1Internal(Long[][] G, double tau,
			int numOfIterations, long m_A, long m_B, Boolean isLogged) throws SFTException{
		
		// adapt parameters to direct product call
		Long[] directG = SFTUtils.getDirectProdGFromAbelianG(G);
		
		// call the corresponding direct product method
		MatlabTemporaryRepositoryDirectProd rep = runMatlabSFTPart1Internal(isLogged,directG,tau,m_A,m_B);
		
		// create a finite-abelian corresponding repository and return it
		return SFTUtils.getMatlabFiniteAbelianRep(rep, G);
	}
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static MatlabTemporaryRepositoryFiniteAbelian runMatlabSFTPart1Internal(Long[][] G, double delta_t, double tau,
			int numOfIterations, double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff, Boolean isLogged)
	throws SFTException{
		
		// adapt parameters to direct product call
		Long[] directG = SFTUtils.getDirectProdGFromAbelianG(G);
		
		// call the corresponding direct product method
		MatlabTemporaryRepositoryDirectProd rep = runMatlabSFTPart1Internal(isLogged,directG,tau,
				delta_t,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);
		
		// create a finite-abelian corresponding repository and return it
		return SFTUtils.getMatlabFiniteAbelianRep(rep, G);
	}
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static MatlabTemporaryResultFiniteAbelian runMatlabSFTPart2Internal(Long[][] G, double tau,
			MatlabTemporaryRepositoryFiniteAbelian matlabRep,int numOfIterations) throws SFTException,FunctionException{
		
		// adapt parameters to direct product call
		Long[] directG = SFTUtils.getDirectProdGFromAbelianG(G);
		MatlabTemporaryRepositoryDirectProd directRep = SFTUtils.getMatlabDirectProdRep(matlabRep, G);
		
		// call the corresponding direct product method, translate into finite Abelian representation and return
		return SFTUtils.getMatlabFiniteAbelianRes(
				runMatlabSFTPart2Internal(directG,tau,directRep,numOfIterations),G);

	}
	
	/* #########################################################################################
	 * 
	 * 							Implementation of the SFT algorithm
	 * 										For JAVA Usage
	 * 
	 * #########################################################################################*/
	
	/* ***********************************************************************
	 * The SFT algorithms 3.9 - 3.12 implementation (as described in the paper)
	 *************************************************************************/
	/**
	 * Main SFT procedure (3.9)
	 * The main SFT is departed into two parts, where part one builds a set of elements to be
	 * f-valued, and part two continues its calculations using these query results.
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned.
	 * @param func
	 * 				The given function over G -> C whose elements and their significant coefficients are returned.
	 * 				Used for query access.
	 * @param numOfIterations
	 * 				The number of SFT procedure iterations to run. Each iteration is ran with the difference function
	 * 				of the given function and the output of the previous SFT iteration.
	 * 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
	 * 				with greater precision.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 */
	@SuppressWarnings("unchecked")
	protected static Map<long[],Complex> runMainSFTAlgorithm(long[] G, double tau, DirectProdFunction func,
			int numOfIterations, long m_A, long m_B) throws SFTException,FunctionException{
		Log.log("SFT -> runMainSFTAlgorithm - main algorithm started");
		
		// call part 1
		Set<long[]>[][] sets = new HashSet[G.length+1][];
		Set<long[]> Q = new HashSet<long[]>();
		callPart1(G,tau,m_A,m_B,sets,Q);
		
		// call part 2
		Map<long[],Complex> res = new HashMap<long[],Complex>();
		callPart2(G,tau,func,numOfIterations,sets,Q,res);
		
		Log.log("SFT -> runMainSFTAlgorithm  - main algorithm completed");
		return res;
	}
	
	/* *****************************************************
	 * implementation of algorithm 3.9 departed into 2 parts
	 * for usage both in Java and in Matlab
	 * *****************************************************/
	/**
	 * Algorithm 3.9 part 1
	 */
	private static void callPart1(long[] G, double tau, long m_A, long m_B,
			Set<long[]>[][] sets, Set<long[]> Q){
		generateQueries(G, m_A, m_B, sets);
		Log.log("\tgenerated sets A,B1,..,BNt for t in {1,...,k} ");
		
		// Build set Q
		//Q = new HashSet<long[]>();
		Set<long[]> A = sets[0][0];
		
		int qSize = 0;
		int k = G.length;
		for (int t=1; t<=k; t++){
			for (int l=0; l<sets[t].length; l++){
				qSize += A.size() * sets[t][l].size();
			}
		}
		Log.log("\tstarting to generate set Q of maximum size of "+qSize);
		
		long qCalcCounter = 0;
		for (int t=1; t<=k; t++){
			Set<long[]>[] tSets = sets[t];
			for (int l=0; l<tSets.length; l++){
				Set<long[]> Btl = tSets[l];
				for(long[] e_a: A){
					for(long[] e_b: Btl){
						long[] elem = SFTUtils.subVectorModulo(e_a, e_b, G[t-1],k); // vector subtraction modulo Nt
						if (!SFTUtils.contains(Q,elem))
							Q.add(elem);
						qCalcCounter++;
						if (qCalcCounter % 10000 == 0)
							Log.log("\tCalculating Q, already checked "+qCalcCounter+" couples of a in A, b in Btl");
					}
				}
			}
		}
		Log.log("\tdone calculating Q, actual size is "+Q.size());
				
		if (Q.size() < 3000){
			String QValues = "";
			int rowCount = 0;
			for (Iterator<long[]> j = Q.iterator(); j.hasNext();){
				rowCount++;
				QValues += SFTUtils.vectorToString(j.next());
				QValues += (rowCount % 20 == 0)? "\n\t":" ";
			}
			String Qsize = Q.size()+"";
			Log.log("\tcreated Q = {x - y | x in A, y in union(B_t_l), t=1,...,k, l=1,...,log(Nt)} of size "+Qsize+":\n\t"+QValues);
		} else Log.log("\tQ is to big to print...");
	}
	
	/**
	 * Algorithm 3.9 part 2
	 */
	private static void callPart2(long[] G,double tau,DirectProdFunction func,int numOfIterations,Set<long[]>[][] sets,
			Set<long[]> Q,Map<long[],Complex> res) throws FunctionException{
		/*
		 * Run iterations of the following:
		 * - estimate the function of the current iteration f_i
		 * - calculate f_(i+1) = f - f_i
		 */
		
		// define result map, accumulates elements and their coefficients along the iterations
		//res = new HashMap<long[],Complex>();
		// allow function access to the accumulating result
		SFTUtils.ResultFunction resFunc = new SFTUtils.ResultFunction(G, res);
		// allow access to the difference function [ f - resFunc ]
		// initial value [ f - 0 ] = f
		SFTUtils.DiffFunction diffFunc = new SFTUtils.DiffFunction(G, func, resFunc);
		
		int i;
		for (i=1; i<=numOfIterations; i++){
			Log.log("\t--- Starting iteration "+i+" out of "+numOfIterations+" ---");
			// create query for current iteration
			Map<String,Complex> query = new HashMap<String,Complex>();
			for(long[] elem: Q){
				query.put(SFTUtils.vectorToString(elem), diffFunc.getValue(elem)); // diffFunc updates from iteration to iteration
			}
			Log.log("\t\tCreated query of size "+query.size());

			// run getFixedQueriesSFT to get the current iteration's elements with significant coefficients
			Set<long[]> tmpRes = getFixedQueriesSFT(G,tau,sets,query);
			Log.log("\t\tcurrent iteration's L size is: "+tmpRes.size());

			// calculate the coefficients for the elements in tmpRes
			Map<long[],Complex> currRes = SFTUtils.calcElemCoeffPairs(tmpRes, query, G);
			Log.log("\t\tCalculated coefficients for current iteration L of iteration "+i);
			
			// add current iteration's results into global result map, which will update diffFunc as well
			for(long[] elem:currRes.keySet()) res.put(elem, currRes.get(elem));
			Log.log("\t\tUpdated difference function");
			
			Log.log("\t--- Done with iteration "+i+" ---");
		}
		
		String LValues;
		Set<long[]> L = res.keySet();
		if (L.size() < 3000){
			LValues = "";
			for (long[] e: L){
				LValues += "\t<"+SFTUtils.vectorToString(e)+","+(res.get(e))+">\n";
			}
		} else LValues = "L is too large to print.";
		Log.log("\tfinished calculating L for "+(i-1)+" iterations, the list of significant Fourier coefficients for f:\n"+LValues+"\n");
		Log.log("\tL is of size "+L.size());
	}
	
	/**
	 * Calculates m_A and m_B as described in algorithm 3.10 (generate queries)
	 */
	protected static long[] getRandomSetsSizes(long[] G, double delta_t, double tau, double fInfNorm, double fEuclideanNorm,
			double deltaCoeff, double randSetsCoeff){
		Log.log("\tcalculating m_A, m_B started");
		
		double gamma = tau/36;
		long sizeOfG = 1;
		int k = G.length;
		for (int i=0; i<k; i++) sizeOfG *= G[i];
		double delta = SFTUtils.calcDelta(delta_t,deltaCoeff,fEuclideanNorm,tau,sizeOfG);
		
		double eta = Math.min(Math.min(gamma, Math.sqrt(gamma)),(gamma/fInfNorm));
		double tmpCoeff = Math.pow(fInfNorm/eta, 2);
		long m_A = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(1.0/delta)));
		long m_B = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(fInfNorm/(delta*gamma))));
		
		Log.log("\tcalculating m_A, m_B finished");
		return new long[]{m_A, m_B};
	}
		
	/**
	 * Generate Queries algorithm (3.10)
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				Sets of elements in G from which the main procedure will
	 * 				create the set of x's to ask their f-value.
	 */
	@SuppressWarnings("unchecked")
	protected static void generateQueries(long[] G, long m_A, long m_B,
			Set<long[]>[][] res){
		Log.log("SFT -> generateQueries started");
		Log.log("\tm_A is: "+m_A+", m_B is: "+m_B);
		
		// generate A,B_1,...,B_Ntl for each t in {1,...,k} and l in {1,...,logN_t}
		int[] logG = SFTUtils.calcLogG(G);
		//Set<long[]>[][] res = new HashSet[G.length+1][];
		
		// generate random subset A partial to G with m_A elements
		res[0] = new HashSet[1];
		res[0][0] = SFTUtils.generateRandomSubsetA(m_A, G);
		
		String AValues = "\t\t";
		long[][] Aarr = new long[res[0][0].size()][]; res[0][0].toArray(Aarr);
		for (int j=0; j<Aarr.length; j++){
			AValues += SFTUtils.vectorToString(Aarr[j])+" ";
			if (j % 20 == 0 && j != 0) AValues+="\n\t\t";
		}
		Log.log("\tA: \n"+AValues+"\n\tEnd of A");

		// generate for each t in {1,...,k} logN_t random subsets B_tl partial to 
		// Z_N1 X ... X Z_Nt-1 X {0,...,2^(l-1)-1} X {0} X ... X {0}
		// and if k=1, then partial to {0,...,2^(l-1)-1} with min{m_B,2^(l-1)} elements
		// return an array of A and B1,...,BlogN_t for each t in {1,...,k}
		for (int t=1; t<=G.length; t++) res[t] = new HashSet[logG[t-1]];
		for (int t=1; t<=G.length; t++)
			for(int l=1; l<=logG[t-1]; l++){
				res[t][l-1] = SFTUtils.generateRandomSubsetBtl(m_B,G,t,l);
			}
		
		Log.log("\tB's:\n");
		int t,l;
		for (t=1; t<=G.length; t++){
			for(l=0; l<logG[t-1]; l++){
				String BtlValues = "size: "+res[t][l].size()+"; elements:\n\t\t";
				long[][] Btlarr = new long[res[t][l].size()][]; res[t][l].toArray(Btlarr);
				for (int j=0; j<Btlarr.length; j++){
					BtlValues += SFTUtils.vectorToString(Btlarr[j])+" ";
					if (j % 20 == 0 && j != 0) BtlValues+="\n\t\t";
				}
				Log.log("\n\tB_"+t+"_"+(l+1)+": "+BtlValues);
			}
		}
		Log.log("\tEnd of B's");
		
		Log.log("\tcreated A and and B1,...,BlogN_t for each t in {1,...,k}");
		Log.log("SFT -> generateQueries completed");
	}
	
	/**
	 * Fixed Queries SFT algorithm (3.11)
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned.
	 * @param querySets
	 * 				The output of the generateQueries function.
	 * @param query
	 * 				The mapping {q,f(q)}
	 * @return
	 * 				A list L of vectors in G of the tau-significant Fourier coefficients
	 * 				of f with probability as set by the initial parameters.
	 */
	protected static Set<long[]> getFixedQueriesSFT(long[] G, double tau, Set<long[]>[][] querySets, Map<String,Complex> query){
		Log.log("SFT -> getFixedQueriesSFT started");
		
		FileWriter f = null; PrintWriter p = null; // TODO remove this shit
		/*
		try{
			f = new FileWriter("matlab\\wav\\est_dist.txt");
			p = new PrintWriter(f);
		} catch (IOException ioe){
			System.err.println("You SUCK!!! IO thrown on writing est_dist.txt");
		}
		*/
		
		int k = G.length;
		Set<long[]> A = querySets[0][0];
		Set<long[]> prefixes = new HashSet<long[]>();
		long[] firstPrefixVec = new long[]{-1}; // first vector for running the first iteration once
		prefixes.add(firstPrefixVec);
		
		// run iterations over t = 1,...,k
		for(int t=1; t<=k; t++){
			
			Set<long[]> tmpPrefixes = new HashSet<long[]>();
			long N = G[t-1];
			// run iterations over l = 0,...,log_2(N)-1
			int logN = SFTUtils.calcLogN(N);
			
			Log.log("\t\t>>> Prefix vectors for stage t = "+t+" (total of "+prefixes.size()+" vectors):");
			String prefixVectorsString = "";
			for (long[] prefixVector: prefixes) prefixVectorsString += SFTUtils.vectorToString(prefixVector)+" ";
			Log.log("\t\t>>> "+prefixVectorsString);
			
			int vecCount = 1;
			for (long[] prefixVector: prefixes){
				Log.log("\t\t\t"+(vecCount++)+" of "+prefixes.size()+" - doing calculation for prefix vector "+SFTUtils.vectorToString(prefixVector)+"...");
				// check if this is the first iteration, that is t == 1. if so, mark the prefix vector as
				// "the empty string", which will be null.
				if (t == 1) prefixVector = null;
				
				// initialize candidate (candidate_0)
				long[] initInterval = {0,N};
				Candidate candidate = new Candidate(initInterval);

				for(int l=0; l<logN; l++){
					Candidate tmpCandidate = new Candidate();
					for (long[] interval: candidate.getSet()){
						// create two sub intervals
						long a = interval[0];
						long b = interval[1];

						long middle = (long)Math.floor((a+b)/2);
						long[] subInterval1 = {a, middle};
						long[] subInterval2 = {middle+1, b};

						// check that the intervals size difference is no greater that 1
						assert(Math.abs((subInterval1[1]-subInterval1[0]) -	(subInterval2[1]-subInterval2[0])) <= 1);
						
						// for each sub interval check if it is "heavy"
						Set<long[]> B_t_lplus1 = querySets[t][l];
						//String vec = (prefixVector == null) ? "empty string" : SFTUtils.vectorToString(prefixVector); 
						//Debug.log("\tcalling distinguish for prefix "+vec+", t="+t+", l="+(l+1)+":");
						//Debug.log("\tsub-interval ["+a+","+middle+"]:");
						if (distinguish(prefixVector, k, G, N, subInterval1, tau, A, B_t_lplus1, query, p))
							tmpCandidate.addInterval(subInterval1);
						//Debug.log("\tsub-interval ["+middle+","+b+"]:");
						if (distinguish(prefixVector, k, G, N, subInterval2, tau, A, B_t_lplus1, query, p))
							tmpCandidate.addInterval(subInterval2);
					}
					candidate = tmpCandidate; // update candidate_i to candidate_(i+1)
				}
				
				// build Lt_prefix for this prefixVector and add it to the prefixes set of the current t
				Set<long[]> Lt_prefix = new HashSet<long[]>();
				// first time
				if (prefixVector == null){
					for (long[] interval: candidate.getSet()){
						if (interval[0] == interval[1]){
							long[] elem = new long[]{interval[0]};
							Lt_prefix.add(elem);
						}
					}
				}
				// not first time
				else {
					for (long[] interval: candidate.getSet()){
						if (interval[0] == interval[1]){
							long[] elem = new long[t];
							// fill first t-1 places with the prefix and last place with the candidate
							int i = 0;
							for (;i<t-1;i++) elem[i] = prefixVector[i];
							elem[i] = interval[0];
							Lt_prefix.add(elem);
						}
					}
				}
				// add all to t's prefix
				tmpPrefixes.addAll(Lt_prefix);
			}
			
			// update prefixes and move to the next t
			prefixes = tmpPrefixes;
		}
		
		Log.log("\tcandidate iterations finished");
		
		Log.log("\tDone creating L");
		Log.log("SFT -> getFixedQueriesSFT completed");
		
		try{ //TODO remove this shit also
			if (f != null && p != null){p.close(); f.close();};
		} catch (IOException ioe){
			System.err.println("closing est_dist.txt failed");
		}
		return prefixes;
	}
	
	/**
	 * Distinguishing algorithm (3.12)
	 * @param prefixVector
	 * 			current prefix vector for the distinguish calculation.
	 * 			If this is the first iteration, that is t == 1, will be null.
	 * @param k
	 * 			the size of G
	 * @param N
	 * 			the size of Z_Nt for the current t
	 * @param interval
	 * 			the interval to be checked for "heaviness"
	 * @param tau
	 * 			threshold
	 * @param A		
	 * @param B
	 * @param query
	 * @return
	 * 			Decision whether to keep or discard the interval {a,b} 
	 */
	protected static boolean distinguish(long[] prefixVector, int k, long[] G, long N, long[] interval, double tau,
			Set<long[]> A, Set<long[]> B, Map<String,Complex> query, PrintWriter p){
		//Debug.log("SFT -> distinguish started");
		
		double est = 0;
		long v = (long)-Math.floor((interval[0]+interval[1])/2);
		int tIndex = (prefixVector == null)?0:prefixVector.length;
		
		// calculate est(a,b)
		for (long[] x: A){
			Complex tmpBSum = new Complex(0,0);
			for (long[] y: B){
				String x_sub_y = SFTUtils.vectorToString(SFTUtils.subVectorModulo(x, y, N, k));
				Complex chiValue = SFTUtils.chi(N,v,y[tIndex]);
				// make it conjucate
				chiValue = chiValue.getConjugate();
				if (prefixVector != null){
					Complex prefixChiValue = SFTUtils.chi(tIndex,G,prefixVector,y);
					chiValue = Complex.mulComplex(chiValue, prefixChiValue);
				}
				tmpBSum.addComplex(Complex.divComplex(Complex.mulComplex(chiValue,query.get(x_sub_y)),(double)B.size()));
			}
			
			est += tmpBSum.getNormSquare()/(double)A.size();
		}
		
		// compare to threshold and return result
		double threshold = 5*tau/36;
		Log.log("\tcalculated est:"+est+((est >= threshold) ? "\t\tPASSED!":""));
		if (interval[0] == interval[1] && p != null) p.println(interval[0]+" "+est);
		
		//Log.log("SFT -> distinguish completed");
		
		return est >= threshold;
	}
	
	/* #########################################################################################
	 * 
	 * 							Implementation of the SFT algorithm
	 * 									 For MATLAB Usage
	 * 
	 * #########################################################################################*/
	
	/* *************************************************************
	 * The SFT algorithm 3.9 departed into 2 parts for use in Matlab
	 ***************************************************************/
	
	/**
	 * Main SFT procedure (3.9) - Part 1/2 (for use in Matlab)
	 * The main SFT is departed into two parts, where part one builds a set of elements to be
	 * f-valued, and part two continues its calculations using these query results.
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned.
	 * @param numOfIterations
	 * 				The number of SFT procedure iterations to run. Each iteration is ran with the difference function
	 * 				of the given function and the output of the previous SFT iteration.
	 * 				This is an optimization for the original SFT algorithm to enable catching significant coefficients
	 * 				with greater precision.
	 * @param m_A
	 * 				The size of the group A (for constructing the group Q).
	 * @param m_B
	 * 				The size of the groups Btl (for constructing the group Q).
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 */
	@SuppressWarnings("unchecked")
	private static Set<long[]>[][] runMatlabSFTPart1Internal(long[] G, double tau, long m_A, long m_B)
	throws SFTException{
		Log.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 1 started");

		// call part 1
		Set<long[]>[][] sets = new HashSet[G.length+1][];
		Set<long[]> Q = new HashSet<long[]>();
		callPart1(G,tau,m_A,m_B,sets,Q);
		
		// put Q as the last set in sets and return A, all B's and Q
		Set<long[]>[][] res = new HashSet[sets.length+1][];
		int i;
		for(i=0; i<sets.length; i++){
			res[i] = new HashSet[sets[i].length];
			for(int j=0; j<sets[i].length; j++){
				res[i][j] = sets[i][j];
			}
		}
		res[i] = new HashSet[1];
		res[i][0] = Q;
		return res;
	}
	
	/**
	 * Main SFT procedure (3.9) - Part 2/2 (for use in Matlab)
	 * The main SFT is departed into two parts, where part one builds a set of elements to be
	 * f-valued, and part two continues its calculations using these query results.
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned.
	 * @param sets
	 * 				The sets A,Btl,Q.
	 * @param query
	 * 				The mapping {q,f(q)}.
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				confidence set by the selection of the m_A and m_B values.
	 */
	private static Map<long[],Complex> runMatlabSFTPart2Internal(long[] G, double tau, Set<long[]>[][] sets,
			Map<String,Complex> query_in, int numOfIterations) throws SFTException,FunctionException{	
		Log.log("SFT -> runMatlabSFTPart2Internal - main algorithm part 2 started");
		
		// create func
		DirectProdFunction func = new SFTUtils.ResultFunction(G, query_in, 1);
		Set<long[]> Q = sets[sets.length-1][0];
		
		// call part 2
		Map<long[],Complex> res = new HashMap<long[],Complex>();
		callPart2(G,tau,func,numOfIterations,sets,Q,res);
		
		Log.log("SFT -> runMainSFTAlgorithm  - main algorithm completed");
		return res;
	}
}
