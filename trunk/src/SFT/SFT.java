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


/**
 * @author Elizabeth Firman and Ariel Stolerman
 * The SFT class provides static methods for calculating the elements in Z_N whose coefficients in a given
 * function over Z_N -> C are significant, done in a 
 * 
 * The implementation is based on the SFT algorithm described in "Learning Noisy Characters, Multiplication Codes, 
 * And Cryptographic Hardcore Predicates" (Adi Akavia, 2008, page 52).
 * 
 * This library is a project in CS workshop, TAU, Spring 2010.
 */
public class SFT {
	
	/* **************************
	 * Interface public functions
	 ****************************/
	
	/**
	 * Returns a set of the elements in Z_N whose coefficients in the given function are tau-significant with
	 * delta-confidence.
	 * @param N
	 * 				The value describing the group Z_N.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over Z_N -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access.
	 * @return
	 * 				A set of the elements in Z_N whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				TODO
	 */
	public static Set<Long> getSignificatElements(long N, double delta, double tau, Function func)
	throws SFTException{
		return null;
	}
	
	/**
	 * Returns a set of the elements in Z_N whose coefficients in the given function are tau-significant with
	 * delta-confidence.
	 * The algorithm includes a calculation of the error-bound, based on the delta-input and a some constant.
	 * This implementation allows the user (who knows the algorithm) to state this constant.
	 * @param N
	 * 				The value describing the group Z_N.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over Z_N -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access. 
	 * @param fInfNorm
	 * 				The infinity norm of the function.
	 * @param fEuclideanNorm
	 * 				The Euclidean norm of the function.
	 * @return
	 * 				A set of the elements in Z_N whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				TODO
	 */
	public static Set<Long> getSignificatElements(long N, double delta, double tau, Function func,
			double fInfNorm, double fEuclideanNorm) throws SFTException{
		return null;
	}
	
	/**
	 * Returns a set of the elements in Z_N whose coefficients in the given function are tau-significant with
	 * delta-confidence.
	 * The algorithm also includes a calculation of log(N)+1 randomly generated sets of elements in Z_N,
	 * of sizes defined as m_A and m_B in the paper, that uses some constant. This implementation allows the user
	 * (who knows the algorithm) to state this constant as well.
	 * @param N
	 * 				The value describing the group Z_N.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over Z_N -> C whose Fourier coefficients (elements) are returned.
	 * 				Used for query access. 
	 * @param deltaCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @param randSetsCoeff
	 * 				A constant coefficient for the algorithm's calculation of delta.
	 * @return
	 * 				A set of the elements in Z_N whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				TODO
	 */
	public static Set<Long> getSignificatElements(long N, double delta, double tau, Function func,
			float deltaCoeff, float randSetsCoeff) throws SFTException{
		return null;
	}	
	
	/**
	 * Returns a set of the elements in Z_N whose coefficients in the given function are tau-significant with
	 * delta-confidence.
	 * The algorithm includes a calculation of the error-bound, based on the delta-input and a some constant.
	 * This implementation allows the user (who knows the algorithm) to state this constant.
	 * The algorithm also includes a calculation of log(N)+1 randomly generated sets of elements in Z_N,
	 * of sizes defined as m_A and m_B in the paper, that uses some constant. This implementation allows the user
	 * (who knows the algorithm) to state this constant as well.
	 * @param N
	 * 				The value describing the group Z_N.
	 * @param delta
	 * 				The confidence parameter such that the algorithm succeeds with probability 1-delta.
	 * @param tau
	 * 				The threshold such that all tau-significant elements are returned. 
	 * @param func
	 * 				The given function over Z_N -> C whose Fourier coefficients (elements) are returned.
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
	 * 				A set of the elements in Z_N whose coefficients in the given function are tau-significant
	 * 				with delta-confidence.
	 * @throws SFTException
	 * 				TODO
	 */
	public static Set<Long> getSignificatElements(long N, double delta, double tau, Function func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		// check inputs (will throw SFTException on input error)
		//SFTUtils.checkParameters(N, delta, tau, fInfNorm, fEuclideanNorm, deltaCoeff, randSetsCoeff);
		// run the SFT algorithm
		//TODO
		return null;
	}
	
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
	protected static Set<long[]> runMainSFTAlgorithm(long[] G, double delta_t, double tau, Function func,
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
		
		// create query
		Map<String,Complex> query = new HashMap<String,Complex>();
		for(long[] elem: Q){
			query.put(SFTUtils.vectorToString(elem), func.getValue(elem));
		}
		
		Debug.log("\tCreated query");
		
		// run getFixedQueriesSFT and return its output, L
		Set<long[]> L = getFixedQueriesSFT(G,tau,sets,query);
		
		String LValues = "";
		for (long[] e: L){
			LValues += SFTUtils.vectorToString(e)+"\n";
		}
		Debug.log("\tfinished calculating L, the list of significant Fourier coefficients for f: "+LValues);
		
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
		
		// generate A,B_1,...,B_Ntl for each t in {1,...,k} and l in {1,...,logN_t}
		
		int[] logG = SFTUtils.calcLogG(G);
		Set<long[]>[][] res = new HashSet[G.length+1][];
		
		// generate random subset A partial to G with m_A elements
		res[0] = new HashSet[1];
		res[0][0] = SFTUtils.generateRandomSubsetA(m_A, G);
		
		String AValues = "";
		for (Iterator<long[]> j = res[0][0].iterator(); j.hasNext(); ){
			AValues += "\t\t"+SFTUtils.vectorToString(j.next())+"\n";
		}
		Debug.log("\tA: \n"+AValues+"End of A");

		// generate for each t in {1,...,k} logN_t random subsets B_tl partial to 
		// Z_N1 X ... X Z_Nt-1 X {0,...,2^(l-1)-1} X {0} X ... X {0}
		// and if k=1, then partial to {0,...,2^(l-1)-1} with min{m_B,2^(l-1)} elements
		// return an array of A and B1,...,BlogN_t for each t in {1,...,k}
		System.out.println(G.length+" "+res.length);
		for (int t=1; t<=G.length; t++) res[t] = new HashSet[logG[t-1]];
		for (int t=1; t<=G.length; t++)
			for(int l=1; l<=logG[t-1]; l++)
				res[t][l-1] = SFTUtils.generateRandomSubsetBtl(m_B,G,t,l);
				
		Debug.log("\tB's:\n");
		int t,l;
		for (t=1; t<=G.length; t++){
			for(l=0; l<logG[t-1]; l++){
				String BtlValues = "size: "+res[t][l].size()+"; elements:\n\t\t";
				for (Iterator<long[]> j = res[t][l].iterator(); j.hasNext(); ){
					BtlValues += SFTUtils.vectorToString(j.next())+" ";
					if (BtlValues.length() % 100 == 0) BtlValues+="\n\t\t";
				}
				Debug.log("\tB_"+t+"_"+(l+1)+": "+BtlValues);
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
			
			for (long[] prefixVector: prefixes){
				// check if this is the first iteration, that is t == 1. if so, mark the prefix vector as
				// "the empty string", which will be null.
				if (t == 1) prefixVector = null;
				
				// initialize candidate (candidate_0)
				long[] initInterval = {0,N}; //TODO: possibly should be N-1, since N is not a candidate!
				Candidate candidate = new Candidate(initInterval);

				// run iterations over l = 0,...,log_2(N)-1
				int logN = SFTUtils.calcLogN(N);

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
						Set<long[]> B_lplus1 = querySets[t][l];
						if (distinguish(prefixVector, k, G, N, subInterval1, tau, A, B_lplus1, query))
							tmpCandidate.addInterval(subInterval1);
						if (distinguish(prefixVector, k, G, N, subInterval2, tau, A, B_lplus1, query))
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
	protected static boolean distinguish(long[] prefixVector, int k, long[] G, long N, long[] interval, double tau, Set<long[]> A,
			Set<long[]> B, Map<String,Complex> query){
		Debug.log("SFT -> distinguish started");
		
		double est = 0;
		long v = (long)-Math.floor((interval[0]+interval[1])/2);
		int tIndex = (prefixVector == null)?0:prefixVector.length-1;
		
		// calculate est(a,b)
		for (long[] x: A){
			double tmpBSum = 0;
			for (long[] y: B){
				String x_sub_y = SFTUtils.vectorToString(SFTUtils.subVectorModulo(x, y, N, k));
				Complex chiValue = SFTUtils.chi(N,v,y[tIndex]);
				if (prefixVector != null)
					chiValue = Complex.mulComplex(chiValue, SFTUtils.chi(tIndex-1,G,prefixVector,y));
				tmpBSum += SFTUtils.innerProduct(chiValue,query.get(x_sub_y));
			}
			tmpBSum /= B.size();
			tmpBSum *= tmpBSum;
			est += tmpBSum;
		}
		
		est /= A.size();
		
		Debug.log("\tcalculated est: "+est);
		
		// compare to threshold and return result
		double threshold = 5*tau/36;
		
		Debug.log("SFT -> distinguish completed");
		
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
		res[i][1] = Q;
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
			LValues += SFTUtils.vectorToString(e)+"\n";
		}
		Debug.log("\tfinished calculating L, the list of significant Fourier coefficients for f: "+LValues);
		
		Debug.log("SFT -> runMatlabSFTPart1Internal - main algorithm part 2 completed");
		return L;
	}
	
}
