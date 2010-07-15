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

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
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
	
	/* *******************************************
	 * delta calculation and random sets constants
	 * *******************************************/
	private static double deltaCoeff = 1;
	private static double randSetsCoeff = 0.0001;
	
	/* *******************************************************************
	 * Interface public functions for direct product G (Z_N1 x ... x Z_Nk)
	 *********************************************************************/
	
	/**
	 * Returns a map of the elements in G and their tau-significant coefficients in the given function with
	 * delta-confidence.
	 * @param G
	 * 				The values N1,...,Nk describing the group G = Z_N1 x ... x Z_Nk.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
	 * @param fInfNorm
	 * 				The infinity norm of the function.
	 * @param fEuclideanNorm
	 * 				The Euclidean norm of the function.
	 * @return
	 * 				A map of the elements in G and their tau-significant coefficients in the given function with
	 * 				delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 */
	public static Map<long[],Complex> getSignificantElements(long[] G, double delta, double tau, DirectProdFunction func,
			double fInfNorm, double fEuclideanNorm,int numOfIterations) throws SFTException,FunctionException,IOException{
		return runMainSFTAlgorithm(G, delta, tau, func, fInfNorm, fEuclideanNorm, deltaCoeff, randSetsCoeff,numOfIterations);
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
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
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
	 */
	public static Map<long[],Complex> getSignificantElements(long[] G, double delta, double tau, DirectProdFunction func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff, int numOfIterations)
			throws SFTException,FunctionException,IOException{
		return runMainSFTAlgorithm(G, delta, tau, func, fInfNorm, fEuclideanNorm, deltaCoeff, randSetsCoeff,numOfIterations);
	}
	
	/* *****************************************************
	 * Interface public functions for finite Abelian G
	 *******************************************************/

	/**
	 * Returns a set of the elements in G whose coefficients in the given function are tau-significant with
	 * delta-confidence.
	 * @param G
	 * 				The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
	 * 				corresponding generators for Nj.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
	 * @param fInfNorm
	 * 				The infinity norm of the function.
	 * @param fEuclideanNorm
	 * 				The Euclidean norm of the function.
	 * @return
	 * 				A set of the elements in G whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 */
	public static Map<Long,Complex> getSignificantElements(long[][] G, double delta, double tau, FiniteAbelianFunction func,
			double fInfNorm, double fEuclideanNorm, int numOfIterations) throws SFTException,IOException{
		// create parameters for the direct product version
		long[] dpG = SFTUtils.getGFromAbelianFunc(func);
		try{
			DirectProdFunction dpFunc = new DirectedProdFromAbelianFunc(dpG,func);
			// call the direct product version of this method
			Map<long[],Complex> Ltag = getSignificantElements(dpG,delta,tau,dpFunc,fInfNorm,fEuclideanNorm,numOfIterations);
			// return the finite Abelian representation of the result set
			return SFTUtils.calcElemCoeffPairsAbelian(Ltag, G);
		} catch (FunctionException fe){
			System.err.println("Invalid function.");
			return null;
		}
	}
	
	/**
	 * Returns a set of the elements in G whose coefficients in the given function are tau-significant with
	 * delta-confidence.
	 * The algorithm includes a calculation of the error-bound, based on the delta-input and a some constant.
	 * This implementation allows the user (who knows the algorithm) to state this constant.
	 * The algorithm also includes a calculation of randomly generated sets of elements in G,
	 * of sizes defined as m_A and m_B in the paper, that uses some constant. This implementation allows the user
	 * (who knows the algorithm) to state this constant as well.
	 * @param G
	 * 				The values (g1,N1),...,(gk,Nk) describing the Abelian group G where gj are the
	 * 				corresponding generators for Nj.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over G -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
	 * @param fInfNorm
	 * 				The infinity norm of the function.
	 * @param fEuclideanNorm
	 * 				The Euclidean norm of the function.
	 * @param deltaCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @param randSetsCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @return
	 * 				A set of the elements in G whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 */
	public static Map<Long,Complex> getSignificantElements(long[][] G, double delta, double tau, FiniteAbelianFunction func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff, int numOfIterations)
			throws SFTException,IOException{
		// create parameters for the direct product version
		long[] dpG = SFTUtils.getGFromAbelianFunc(func);
		try{
			DirectProdFunction dpFunc = new DirectedProdFromAbelianFunc(dpG,func);
			// call the direct product version of this method
			Map<long[],Complex> Ltag = getSignificantElements(dpG,delta,tau,dpFunc,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff,numOfIterations);
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
	public static MatlabTemporaryRepositoryDirectProd runMatlabSFTPart1Internal(Long[] G, double delta_t, double tau,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff, Boolean isLogged) throws SFTException{
		// set variables to fit algorithm
		Log.setLogMode(isLogged);
		long[] g = new long[G.length];
		for(int i=0; i<G.length; i++) g[i] = G[i];
		
		// run algorithm
		Set<long[]>[][] sets = runMatlabSFTPart1Internal(g, delta_t, tau, fInfNorm, fEuclideanNorm, deltaCoeff, randSetsCoeff);
		
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
	public static Long[][] runMatlabSFTPart2Internal(Long[] G, double tau,
			MatlabTemporaryRepositoryDirectProd matlabRep) throws SFTException, IOException{
		// fit parameters to java
		long[] g = new long[G.length];
		for(int i=0; i<G.length; i++) g[i] = G[i];
		Set<long[]>[][] sets = matlabRep.getSets();
		Map<String,Complex> query = (HashMap<String,Complex>)matlabRep.getQuery();

		// run part 2 and return the results fitted to matlab (Long[] elements)
		Set<long[]> res = runMatlabSFTPart2Internal(g, tau, sets, query);
		Long[][] L = new Long[res.size()][];
		int i=0;
		for (long[] elem: res){
			int len = elem.length;
			Long[] Elem = new Long[len];
			for (int k=0; k<len; k++) Elem[k] = new Long(elem[k]);
			L[i++] = Elem;
		}
		return L;
	}
	
	/* ***********************************************************
	 * Matlab interface public functions for finite Abelian groups
	 * ***********************************************************/
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static MatlabTemporaryRepositoryFiniteAbelian runMatlabSFTPart1Internal(Long[][] G, double delta_t, double tau,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff, Boolean isLogged) throws SFTException{
		
		// adapt parameters to direct product call
		Long[] directG = SFTUtils.getDirectProdGFromAbelianG(G);
		
		// call the corresponding direct product method
		MatlabTemporaryRepositoryDirectProd rep = runMatlabSFTPart1Internal(directG,delta_t,tau,
				fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff,isLogged);
		
		// create a finite-abelian corresponding repository and return it
		return SFTUtils.getMatlabFiniteAbelianRep(rep, G);
	}
	
	/**
	 * For inner use in the Matlab SFT scripts.
	 */
	public static Long[] runMatlabSFTPart2Internal(Long[][] G, double tau,
			MatlabTemporaryRepositoryFiniteAbelian matlabRep) throws SFTException, IOException{
		
		// adapt parameters to direct product call
		Long[] directG = SFTUtils.getDirectProdGFromAbelianG(G);
		MatlabTemporaryRepositoryDirectProd directRep = SFTUtils.getMatlabDirectProdRep(matlabRep, G);
		
		// call the corresponding direct product method
		Long[][] L = runMatlabSFTPart2Internal(directG,tau,directRep);
		
		// create a finite-abelian corresponding L and return it
		Long[] res = new Long[L.length];
		for (int i=0; i<res.length; i++){
			res[i] = SFTUtils.calcAbelianProd(L[i], G);
		}
		return res;
	}
	
	/* #########################################################################################
	 * 
	 * 							Implementation of the SFT algorithm
	 * 
	 * #########################################################################################*/
	
	/* ***********************************************************************
	 * The SFT algorithms 3.9 - 3.12 implementation (as described in the paper)
	 *************************************************************************/
	
	/**
	 * Main SFT procedure (3.9)
	 * The main SFT is departed into two parts, where part one builds a set of elements to be
	 * f-valued, and part two continues its calculations using these query results
	 * @param G			an integer array describing the group Z_N1 X ... X Z_Nk
	 * @param tau		threshold on the weight of the Fourier coefficients we seek
	 * @param delta_t	confidence parameter
	 * @return			a short list L in G of the tau-significant elements and their Fourier coefficients
	 * 					of f with probability at least 1-delta_t
	 */
	protected static Map<long[],Complex> runMainSFTAlgorithm(long[] G, double delta_t, double tau, DirectProdFunction func,
			double fInfNorm, double fEuclideanNorm, double deltaCoeff, double randSetsCoeff, int numOfIterations)
			throws SFTException,FunctionException,IOException{
		Log.log("SFT -> runMainSFTAlgorithm - main algorithm started");
		
		/* run generateQueries (algorithm 3.10) on:
		 * G, gamma = tau/36, ||f||_infinity and delta = delta_t/O((||f||_2^2/tau)^1.5*log_2|G|)
		 */
		double gamma = tau/36;
		long sizeOfG = 0;
		int k = G.length;
		for (int i=0; i<k; i++) sizeOfG += G[i];
		double delta = SFTUtils.calcDelta(delta_t,deltaCoeff,fEuclideanNorm,tau,sizeOfG);
		Log.log("\tgamma is: "+gamma+", delta is: "+delta+", fEuclideanNorm is: "+fEuclideanNorm+", fInfNorm is: "+fInfNorm);
		
		Set<long[]>[][] sets = generateQueries(G, gamma, fInfNorm, delta, randSetsCoeff);
		Log.log("\tgenerated sets A,B1,..,BNt for t in {1,...,k} ");
		
		// Build set Q
		Set<long[]> Q = new HashSet<long[]>();
		Set<long[]> A = sets[0][0];
		
		int qSize = 0;
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
		
		/*
		 * Run iterations of the following:
		 * - estimate the function of the current iteration f_i
		 * - calculate f_(i+1) = f - f_i
		 */
		
		Map<long[],Complex> res = new HashMap<long[],Complex>();
		DirectProdFunction currFunc = func;
		int i;
		for (i=1; i<=numOfIterations; i++){
			Log.log("\t--- Starting iteration "+i+" out of "+numOfIterations+" ---");
			// create query
			Map<String,Complex> query = new HashMap<String,Complex>();
			for(long[] elem: Q){
				query.put(SFTUtils.vectorToString(elem), currFunc.getValue(elem));
			}
			Log.log("\t\tCreated query of size "+query.size());

			// run getFixedQueriesSFT to get the current iteration's elements with significant coefficients
			Set<long[]> tmpRes = getFixedQueriesSFT(G,tau,sets,query);
			Log.log("\t\tcurrent iteration's L size is: "+tmpRes.size());

			// calculate the coefficients for the elements in tmpRes
			Map<long[],Complex> currRes = SFTUtils.calcElemCoeffPairs(tmpRes, query, G);
			Log.log("\t\tCalculated coefficients for current iteration L of iteration "+i);
			
			// add current iteration's results into global result map
			for(long[] elem:currRes.keySet()) res.put(elem, currRes.get(elem));
			
			// calculate the difference function
			currFunc = new SFTUtils.DiffFunction(G, func, currFunc);
			Log.log("\t\tCreated difference function");
			
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
		
		Log.log("SFT -> runMainSFTAlgorithm  - main algorithm completed");
		return res;
	}
		
	/**
	 * Generate Queries algorithm (3.10)
	 * @param G:		an integer array describing the group Z_N1 X ... X Z_Nk
	 * @param gamma:	a value in R+
	 * @param fInfNorm:	||f||_(infinity)
	 * @param deltha:	confidence parameter
	 * @return:			sets of elements in G from which the main procedure will
	 * 					create the set of x's to ask their f-value
	 */
	@SuppressWarnings("unchecked")
	protected static Set<long[]>[][] generateQueries(long[] G, double gamma, double fInfNorm, double delta, double randSetsCoeff){
		Log.log("SFT -> generateQueries started");
		
		// compute m_A and m_B
		double eta = Math.min(Math.min(gamma, Math.sqrt(gamma)),(gamma/fInfNorm));
		double tmpCoeff = Math.pow(fInfNorm/eta, 2);
		long m_A = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(1.0/delta)));
		long m_B = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(fInfNorm/(delta*gamma))));
		m_A = 2*SFTUtils.calcLogN(G[0]);
		m_B = m_A;
		
		Log.log("\tm_A is: "+m_A+", m_B is: "+m_B);
		if (m_A <= 0 || m_B <= 0) System.exit(0);
		//System.exit(0);
		
		// generate A,B_1,...,B_Ntl for each t in {1,...,k} and l in {1,...,logN_t}
		
		int[] logG = SFTUtils.calcLogG(G);
		Set<long[]>[][] res = new HashSet[G.length+1][];
		
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
		
		return res;
	}
	
	/**
	 * Fixed Queries SFT algorithm (3.11)
	 * @param G:			a vector of values describing the G
	 * @param tau:			threshold on the weight of the Fourier coefficients we seek
	 * @param querySets:	the output of the generateQueries function
	 * @param query:		a set {q,f(q)}
	 * @return:				a short list L of vectors in G of the tau-significant Fourier coefficients
	 * 						of f with probability at least 1-deltha_t
	 */
	protected static Set<long[]> getFixedQueriesSFT(long[] G, double tau, Set<long[]>[][] querySets, Map<String,Complex> query) throws IOException{
		Log.log("SFT -> getFixedQueriesSFT started");
		
		FileWriter f = new FileWriter("matlab\\wav\\est_dist.txt");
		PrintWriter p = new PrintWriter(f);
		
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
		
		p.close(); f.close();
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
	 * @return			decides whether to keep or discard the interval {a,b} 
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
		if (interval[0] == interval[1]) p.println(interval[0]+" "+est);
		
		//Log.log("SFT -> distinguish completed");
		
		return est >= threshold;
	}
	
	/* *************************************************************
	 * The SFT algorithm 3.9 departed into 2 parts for use in Matlab
	 ***************************************************************/
	
	/**
	 * Main SFT procedure (3.9) - Part 1/2 (for use in Matlab)
	 * The main SFT is departed into two parts, where part one builds a set of elements to be
	 * f-valued, and part two continues its calculations using these query results.
	 * @param G			an integer array describing the group Z_N1 X ... X Z_Nk
	 * @param tau		threshold on the weight of the Fourier coefficients we seek
	 * @param delta_t	confidence parameter
	 * @return			a short list L in G of the tau-significant Fourier coefficients
	 * 					of f with probability at least 1-delta_t
	 */
	@SuppressWarnings("unchecked")
	private static Set<long[]>[][] runMatlabSFTPart1Internal(long[] G, double delta_t, double tau,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		Log.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 1 started");
		/* run generateQueries (algorithm 3.10) on:
		 * G, gamma = tau/36, ||f||_infinity and delta = delta_t/O((||f||_2^2/tau)^1.5*log_2|G|)
		 */
		double gamma = tau/36;
		long sizeOfG = 0;
		int k = G.length;
		for (int i=0; i<k; i++) sizeOfG += G[i];
		double delta = SFTUtils.calcDelta(delta_t,deltaCoeff,fEuclideanNorm,tau,sizeOfG);
		Log.log("\tgamma is: "+gamma+", delta is: "+delta+", fInfNorm is: "+fInfNorm);
		
		Set<long[]>[][] sets = generateQueries(G, gamma, fInfNorm, delta, randSetsCoeff);
		Log.log("\tgenerated sets A,B1,..,BNt for t in {1,...,k} ");
		
		// Build set Q
		Set<long[]> Q = new HashSet<long[]>();
		Set<long[]> A = sets[0][0];
		
		int qSize = 0;
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
		
		Log.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 1 completed");
		
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
	 * @param G			an integer array describing the group Z_N1 X ... X Z_Nk
	 * @param tau		threshold on the weight of the Fourier coefficients we seek
	 * @param delta_t	confidence parameter
	 * @return			a short list L in G of the tau-significant Fourier coefficients
	 * 					of f with probability at least 1-delta_t
	 */
	private static Set<long[]> runMatlabSFTPart2Internal(long[] G, double tau, Set<long[]>[][] sets,
			Map<String,Complex> query) throws SFTException, IOException{
		Log.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 2 started");
		
		// run getFixedQueriesSFT and return its output, L
		Set<long[]> L = getFixedQueriesSFT(G,tau,sets,query);
		
		String LValues = "";
		for (long[] e: L){
			LValues += "\t"+SFTUtils.vectorToString(e)+"\n";
		}
		Log.log("\tfinished calculating L, the list of significant Fourier coefficients for f:\n"+LValues+"\n");
		
		Log.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 2 completed");
		return L;
	}

	public static double getDeltaCoeff() {
		return deltaCoeff;
	}

	public static double getRandSetsCoeff() {
		return randSetsCoeff;
	}
}
