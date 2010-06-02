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

import java.util.*;
import Function.*;
import SFT.SFTUtils.DirectedProdFromAbelianFunc;
import SFT.SFTUtils.MatlabTemporaryRepository;


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
	
	/* *******************************************************************
	 * Interface public functions for direct product G (Z_N1 x ... x Z_Nk)
	 *********************************************************************/
	
	/**
	 * Returns a set of the elements in G whose coefficients in the given function are tau-significant with
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
	 * 				A set of the elements in G whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 */
	public static Set<long[]> getSignificatElements(long[] G, double delta, double tau, DirectProdFunction func,
			double fInfNorm, double fEuclideanNorm) throws SFTException{
		//TODO
		return null;
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
	 * 				A set of the elements in G whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				If the given parameters are invalid.
	 */
	public static Set<long[]> getSignificatElements(long[] G, double delta, double tau, DirectProdFunction func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		//TODO
		return null;
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
	public static Set<Long> getSignificatElements(long[][] G, double delta, double tau, FiniteAbelianFunction func,
			double fInfNorm, double fEuclideanNorm) throws SFTException{
		// create parameters for the direct product version
		long[] dpG = SFTUtils.getGFromAbelianFunc(func);
		try{
			DirectProdFunction dpFunc = new DirectedProdFromAbelianFunc(dpG,func);
			// call the direct product version of this method
			Set<long[]> Ltag = getSignificatElements(dpG,delta,tau,dpFunc,fInfNorm,fEuclideanNorm);
			// return the finite Abelian representation of the result set
			return SFTUtils.getAbelianRepresentation(Ltag, G);
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
	public static Set<Long> getSignificatElements(long[][] G, double delta, double tau, FiniteAbelianFunction func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		// create parameters for the direct product version
		long[] dpG = SFTUtils.getGFromAbelianFunc(func);
		try{
			DirectProdFunction dpFunc = new DirectedProdFromAbelianFunc(dpG,func);
			// call the direct product version of this method
			Set<long[]> Ltag = getSignificatElements(dpG,delta,tau,dpFunc,fInfNorm,fEuclideanNorm,deltaCoeff,randSetsCoeff);
			// return the finite Abelian representation of the result set
			return SFTUtils.getAbelianRepresentation(Ltag, G);
		} catch (FunctionException fe){
			System.err.println("Invalid function.");
			return null;
		}
	}	
	
	/* *********************************
	 * Matlab interface public functions
	 * ********************************/
	
	public static MatlabTemporaryRepository runMatlabSFTPart1Internal(Long[] G, double delta_t, double tau,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff, Boolean isLogged) throws SFTException{
		// set variables to fit algorithm
		Debug.DEBUG_MODE = isLogged;
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
		MatlabTemporaryRepository matlabRep = new MatlabTemporaryRepository(sets, Q, null);
		return matlabRep;
	}
	
	public static Long[][] runMatlabSFTPart2Internal(Long[] G, double tau,
			MatlabTemporaryRepository matlabRep) throws SFTException{
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
	
	/* ***********************************************************************
	/* ***********************************************************************
	/* ***********************************************************************
	/* ***********************************************************************/
	
	
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
	 * @return			a short list L in G of the tau-significant Fourier coefficients
	 * 					of f with probability at least 1-delta_t
	 */
	protected static Set<long[]> runMainSFTAlgorithm(long[] G, double delta_t, double tau, DirectProdFunction func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		Debug.log("SFT -> runMainSFTAlgorithm - main algorithm started");
		
		/* run generateQueries (algorithm 3.10) on:
		 * G, gamma = tau/36, ||f||_infinity and delta = delta_t/O((||f||_2^2/tau)^1.5*log_2|G|)
		 */
		double gamma = tau/36;
		long sizeOfG = 0;
		int k = G.length;
		for (int i=0; i<k; i++) sizeOfG += G[i];
		double delta = SFTUtils.calcDelta(delta_t,deltaCoeff,fEuclideanNorm,tau,sizeOfG);
		Debug.log("\tgamma is: "+gamma+", delta is: "+delta+", fEuclideanNorm is: "+fEuclideanNorm+", fInfNorm is: "+fInfNorm);
		
		Set<long[]>[][] sets = generateQueries(G, gamma, fInfNorm, delta, randSetsCoeff);
		Debug.log("\tgenerated sets A,B1,..,BNt for t in {1,...,k} ");
		
		// Build set Q
		Set<long[]> Q = new HashSet<long[]>();
		Set<long[]> A = sets[0][0];
		
		int qSize = 0;
		for (int t=1; t<=k; t++){
			for (int l=0; l<sets[t].length; l++){
				qSize += A.size() * sets[t][l].size();
			}
		}
		Debug.log("\tstarting to generate set Q of maximum size of "+qSize);
		
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
							Debug.log("\tCalculating Q, already checked "+qCalcCounter+" couples of a in A, b in Btl");
					}
				}
			}
		}
		Debug.log("\tdone calculating Q, actual size is "+Q.size());
		
		if (Q.size() < 3000){
			String QValues = "";
			int rowCount = 0;
			for (Iterator<long[]> j = Q.iterator(); j.hasNext();){
				rowCount++;
				QValues += SFTUtils.vectorToString(j.next());
				QValues += (rowCount % 20 == 0)? "\n\t":" ";
			}
			String Qsize = Q.size()+"";
			Debug.log("\tcreated Q = {x - y | x in A, y in union(B_t_l), t=1,...,k, l=1,...,log(Nt)} of size "+Qsize+":\n\t"+QValues);
		} else Debug.log("\tQ is to big to print...");
		
		// create query
		Map<String,Complex> query = new HashMap<String,Complex>();
		for(long[] elem: Q){
			query.put(SFTUtils.vectorToString(elem), func.getValue(elem));
		}
		
		Debug.log("\tCreated query of size "+query.size());
		
		// run getFixedQueriesSFT and return its output, L
		Set<long[]> L = getFixedQueriesSFT(G,tau,sets,query);
		Debug.log("\tL size is: "+L.size());
		
		String LValues = "\n";
		for (long[] e: L){
			LValues += "\t"+SFTUtils.vectorToString(e)+"\n";
		}
		Debug.log("\tfinished calculating L, the list of significant Fourier coefficients for f:"+LValues+"\n");
		
		Debug.log("SFT -> runMainSFTAlgorithm  - main algorithm completed");
		return L;
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
	protected static Set<long[]>[][] generateQueries(long[] G, double gamma, double fInfNorm, double delta, float randSetsCoeff){
		Debug.log("SFT -> generateQueries started");
		
		// compute m_A and m_B
		double eta = Math.min(Math.min(gamma, Math.sqrt(gamma)),(gamma/fInfNorm));
		double tmpCoeff = Math.pow(fInfNorm/eta, 2);
		long m_A = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(1.0/delta)));
		long m_B = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(fInfNorm/(delta*gamma))));
		
		Debug.log("\tm_A is: "+m_A+", m_B is: "+m_B);
		if (m_A <= 0 || m_B <= 0) System.exit(0);
		
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
		Debug.log("\tA: \n"+AValues+"\n\tEnd of A");

		// generate for each t in {1,...,k} logN_t random subsets B_tl partial to 
		// Z_N1 X ... X Z_Nt-1 X {0,...,2^(l-1)-1} X {0} X ... X {0}
		// and if k=1, then partial to {0,...,2^(l-1)-1} with min{m_B,2^(l-1)} elements
		// return an array of A and B1,...,BlogN_t for each t in {1,...,k}
		for (int t=1; t<=G.length; t++) res[t] = new HashSet[logG[t-1]];
		for (int t=1; t<=G.length; t++)
			for(int l=1; l<=logG[t-1]; l++){
				res[t][l-1] = SFTUtils.generateRandomSubsetBtl(m_B,G,t,l);
			}
		
		Debug.log("\tB's:\n");
		int t,l;
		for (t=1; t<=G.length; t++){
			for(l=0; l<logG[t-1]; l++){
				String BtlValues = "size: "+res[t][l].size()+"; elements:\n\t\t";
				long[][] Btlarr = new long[res[t][l].size()][]; res[t][l].toArray(Btlarr);
				for (int j=0; j<Btlarr.length; j++){
					BtlValues += SFTUtils.vectorToString(Btlarr[j])+" ";
					if (j % 20 == 0 && j != 0) BtlValues+="\n\t\t";
				}
				Debug.log("\n\tB_"+t+"_"+(l+1)+": "+BtlValues);
			}
		}
		Debug.log("\tEnd of B's");
		
		Debug.log("\tcreated A and and B1,...,BlogN_t for each t in {1,...,k}");
		Debug.log("SFT -> generateQueries completed");
		
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
	protected static Set<long[]> getFixedQueriesSFT(long[] G, double tau, Set<long[]>[][] querySets, Map<String,Complex> query){
		Debug.log("SFT -> getFixedQueriesSFT started");
		
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
			
			Debug.log("\t\t>>> Prefix vectors for stage t = "+t+" (total of "+prefixes.size()+" vectors):");
			String prefixVectorsString = "";
			for (long[] prefixVector: prefixes) prefixVectorsString += SFTUtils.vectorToString(prefixVector)+" ";
			Debug.log("\t\t>>> "+prefixVectorsString);
			
			int vecCount = 1;
			for (long[] prefixVector: prefixes){
				Debug.log("\t\t\t"+(vecCount++)+" of "+prefixes.size()+" - doing calculation for prefix vector "+SFTUtils.vectorToString(prefixVector)+"...");
				// check if this is the first iteration, that is t == 1. if so, mark the prefix vector as
				// "the empty string", which will be null.
				if (t == 1) prefixVector = null;
				
				// initialize candidate (candidate_0)
				long[] initInterval = {0,N}; //TODO: possibly should be N-1, since N is not a candidate!
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
						String vec = (prefixVector == null) ? "empty string" : SFTUtils.vectorToString(prefixVector); 
						//Debug.log("\tcalling distinguish for prefix "+vec+", t="+t+", l="+(l+1)+":");
						//Debug.log("\tsub-interval ["+a+","+middle+"]:");
						if (distinguish(prefixVector, k, G, N, subInterval1, tau, A, B_t_lplus1, query))
							tmpCandidate.addInterval(subInterval1);
						//Debug.log("\tsub-interval ["+middle+","+b+"]:");
						if (distinguish(prefixVector, k, G, N, subInterval2, tau, A, B_t_lplus1, query))
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
		
		Debug.log("\tcandidate iterations finished");
		
		Debug.log("\tDone creating L");
		Debug.log("SFT -> getFixedQueriesSFT completed");
		
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
			Set<long[]> A, Set<long[]> B, Map<String,Complex> query){
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
		//Debug.log("\tcalculated est:"+est+((est >= threshold) ? "\t\tPASSED!":""));
		
		//Debug.log("SFT -> distinguish completed");
		
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
	public static Set<long[]>[][] runMatlabSFTPart1Internal(long[] G, double delta_t, double tau,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		Debug.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 1 started");
		/* run generateQueries (algorithm 3.10) on:
		 * G, gamma = tau/36, ||f||_infinity and delta = delta_t/O((||f||_2^2/tau)^1.5*log_2|G|)
		 */
		double gamma = tau/36;
		long sizeOfG = 0;
		int k = G.length;
		for (int i=0; i<k; i++) sizeOfG += G[i];
		double delta = SFTUtils.calcDelta(delta_t,deltaCoeff,fEuclideanNorm,tau,sizeOfG);
		Debug.log("\tgamma is: "+gamma+", delta is: "+delta+", fInfNorm is: "+fInfNorm);
		
		Set<long[]>[][] sets = generateQueries(G, gamma, fInfNorm, delta, randSetsCoeff);
		Debug.log("\tgenerated sets A,B1,..,BNt for t in {1,...,k} ");
		
		// Build set Q
		Set<long[]> Q = new HashSet<long[]>();
		Set<long[]> A = sets[0][0];
		
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
							Debug.log("\tCalculating Q, already checked "+qCalcCounter+" couples of a in A, b in Btl");
					}
				}
			}
		}
		
		String QValues = "";
		int rowCount = 0;
		for (Iterator<long[]> j = Q.iterator(); j.hasNext();){
			rowCount++;
			QValues += SFTUtils.vectorToString(j.next());
			QValues += (rowCount % 20 == 0)? "\n\t":" ";
		}
		String Qsize = Q.size()+"";
		Debug.log("\tcreated Q = {x - y | x in A, y in union(B_t_l), t=1,...,k, l=1,...,log(Nt)} of size "+Qsize+":\n\t"+QValues);
		
		Debug.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 1 completed");
		
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
	public static Set<long[]> runMatlabSFTPart2Internal(long[] G, double tau, Set<long[]>[][] sets,
			Map<String,Complex> query) throws SFTException{
		Debug.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 2 started");
		
		// run getFixedQueriesSFT and return its output, L
		Set<long[]> L = getFixedQueriesSFT(G,tau,sets,query);
		
		String LValues = "";
		for (long[] e: L){
			LValues += "\n\t"+SFTUtils.vectorToString(e)+"\n";
		}
		Debug.log("\tfinished calculating L, the list of significant Fourier coefficients for f:"+LValues+"\n");
		
		Debug.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 2 completed");
		return L;
	}
}
