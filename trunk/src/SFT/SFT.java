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
		return null;
	}
	
	/* ***********************************************************************
	 * The SFT algorithms 3.4 - 3.7 implementation (as described in the paper)
	 *************************************************************************/
	
	/**
	 * Main SFT procedure (3.4)
	 * The main SFT is departed into two parts, where part one builds a set of elements to be
	 * f-valued, and part two continues its calculations using these query results
	 * @param N:		an integer value describing the group Z_N
	 * @param tau:		threshold on the weight of the Fourier coefficients we seek
	 * @param deltha_t:	confidence parameter
	 * @return			a short list L in Z_N of the tau-significant Fourier coefficients
	 * 					of f with probability at least 1-deltha_t
	 */
	public static void runMainSFTAlgorithm(long N, double delta_t, double tau, Function func,
			double fInfNorm, double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		Debug.log("SFT -> runMainSFTAlgorithm - main algorithm started");
		
		/* run generateQueries (algorithm 3.5) on:
		 * N, gamma = tau/36, ||f||_infinity and delta = delta_t/O((||f||_2^2/tau)^1.5*log_2(N))
		 */
		double gamma = tau/36;
		double delta = SFTUtils.calcDelta(delta_t,deltaCoeff,fEuclideanNorm,tau,N);
		Debug.log("\tgamma is: "+gamma+", delta is: "+delta+", fInfNorm is: "+fInfNorm);
		
		Set<Long>[] sets = generateQueries(N, gamma, fInfNorm, delta, randSetsCoeff);
		Debug.log("\tgenerated sets A,B1,..,Bl");
		
		// Build set Q
		Set<Long> Q = new HashSet<Long>();
		Set<Long> A = sets[0];
		
		for (int i=1; i<sets.length; i++){
			Set<Long> Bl = sets[i];
			for(long e_a: A){
				for(long e_b: Bl){
					long elem = SFTUtils.subModulo(e_a, e_b, N); // subtraction modulo N
					if (!Q.contains(elem)) //TODO check that actually works
						Q.add(elem);
				}
			}
		}
		
		String QValues = "";
		int rowCount = 0;
		for (Iterator<Long> j = Q.iterator(); j.hasNext();){
			rowCount++;
			QValues += j.next();
			QValues += (rowCount % 20 == 0)? "\n\t":" ";
		}
		String Qsize = Q.size()+"";
		Debug.log("\tcreated Q = {x - y | x in A, y in union(B_i), i=1,...,log(N)} of size "+Qsize+":\n\t"+QValues);
		
		// create query
		Map<Long,Complex> query = new HashMap<Long,Complex>();
		for(long elem: Q){
			query.put(elem, func.getValue(elem));
		}
		
		Debug.log("\tCreated query");
		
		// run getFixedQueriesSFT and return its output, L
		Set<Long> L = getFixedQueriesSFT(N,tau,sets,query);
		
		String LValues = "";
		for (long e: L){
			LValues += String.valueOf(e)+" ";
		}
		Debug.log("\tfinished calculating L, the list of significant Fourier coefficients for f: "+LValues);
		
		Debug.log("SFT -> runMainSFTAlgorithmCont  - main algorithm completed");
	}
	
	/**
	 * Generate Queries algorithm (3.5)
	 * @param N:		an integer value describing the group Z_N
	 * @param gamma:	a value in R+
	 * @param fInfNorm:	||f||_(infinity)
	 * @param deltha:	confidence parameter
	 * @return:			sets of elements in Z_N from which the main procedure will
	 * 					create the set of x's to ask their f-value
	 */
	@SuppressWarnings("unchecked")
	public static Set<Long>[] generateQueries(long N, double gamma, double fInfNorm, double delta, float randSetsCoeff){
		Debug.log("SFT -> generateQueries started");
		
		// compute m_A and m_B
		double eta = Math.min(Math.min(gamma, Math.sqrt(gamma)),(gamma/fInfNorm));
		double tmpCoeff = Math.pow(fInfNorm/eta, 2);
		long m_A = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(1.0/delta)));
		long m_B = (long) (randSetsCoeff * Math.ceil(tmpCoeff*Math.log(fInfNorm/(delta*gamma))));
		
		//TODO limit to log(N)?
		//regulation of m_A and m_B - is it correct?
		/*if (m_A > N) {
			Debug.log("\tRegulated m_A from "+m_A+" to "+((long)Math.log(N)));
			m_A = (long)Math.log(N);
		}
		if (m_B > N) {
			Debug.log("\tRegulated m_B from "+m_B+" to "+((long)Math.log(N)));
			m_B = (long)Math.log(N);
		}*/
		
		Debug.log("\tm_A is: "+m_A+", m_B is: "+m_B);
		
		// generate A,B1,...,Bl
		int logN = SFTUtils.calcLogN(N);
		Set<Long>[] res = new HashSet[logN+1];
		
		// generate random subset A partial to Z_N with m_A elements
		res[0] = SFTUtils.generateRandomSubsetA(m_A, N);
		
		String AValues = "";
		for (Iterator<Long> j = res[0].iterator(); j.hasNext(); ){
			AValues += j.next()+" ";
		}
		Debug.log("\tA: "+AValues);

		// generate logN random subsets B_l partial to {0,...,2^(l-1)-1} with min{m_B,2^(l-1)} elements
		// return an array of A,B1,...,Bl
		for(int l=1; l<=logN; l++){
			res[l] = SFTUtils.generateRandomSubsetBl(m_B,N,l);
		}
		
		Debug.log("\tB's:");
		for (int i=1; i<=logN; i++){
			String BlValues = "size: "+res[i].size()+"; elements: ";
			for (Iterator<Long> j = res[i].iterator(); j.hasNext(); ){
				BlValues += j.next()+" ";
			}
			Debug.log("\tB["+i+"]: "+BlValues);
		}
		
		Debug.log("\tcreated A and B1,...,Bl");
		Debug.log("SFT -> generateQueries completed");
		
		return res;
	}
	
	/**
	 * Fixed Queries SFT algorithm (3.6)
	 * @param N:			an integer value describing the group Z_N
	 * @param tau:			threshold on the weight of the Fourier coefficients we seek
	 * @param querySets:	the output of the generateQueries function
	 * @param query:		a set {q,f(q)}
	 * @return:				a short list L in Z_N of the tau-significant Fourier coefficients
	 * 						of f with probability at least 1-deltha_t
	 */
	public static Set<Long> getFixedQueriesSFT(long N, double tau, Set<Long>[] querySets, Map<Long,Complex> query){
		Debug.log("SFT -> getFixedQueriesSFT started");
		
		// initialize candidate (candidate_0)
		long[] initInterval = {0,N}; //TODO: possibly should be N-1, since N is not a candidate!
		Candidate candidate = new Candidate(initInterval);
		
		// run iterations over l = 0,...,log_2(N)-1
		int logN = SFTUtils.calcLogN(N);
		Set<Long> A = querySets[0];
		
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
				Set<Long> B_lplus1 = querySets[l+1];
				if (distinguish(N, subInterval1, tau, A, B_lplus1, query))
					tmpCandidate.addInterval(subInterval1);
				if (distinguish(N, subInterval2, tau, A, B_lplus1, query))
					tmpCandidate.addInterval(subInterval2);
			}
			candidate = tmpCandidate; // update candidate_i to candidate_(i+1)
		}
		
		Debug.log("\tcandidate iterations finished");
		
		// build L and return it
		Set<Long> L = new HashSet<Long>();
		for (long[] interval: candidate.getSet()){
			if (interval[0] == interval[1]){
				long elem = interval[0];
				L.add(elem);
			}
		}
		
		Debug.log("\tDone creating L");
		Debug.log("SFT -> getFixedQueriesSFT completed");
		
		return L;
	}
	
	/**
	 * Distinguishing algorithm (3.7)
	 * @param interval:	the interval to be checked for "heaviness"
	 * @param tau:		threshold
	 * @param A:		
	 * @param B:		
	 * @param query:	
	 * @return:			decides whether to keep or discard the interval {a,b} 
	 */
	public static boolean distinguish(long N, long[] interval, double tau, Set<Long> A, Set<Long> B,
			Map<Long,Complex> query){
		Debug.log("SFT -> distinguish started");
		
		double est = 0;
		long v = (long)-Math.floor((interval[0]+interval[1])/2);
		
		// calculate est(a,b)
		for (long x: A){
			double tmpBSum = 0;
			for (long y: B){
				long x_sub_y = SFTUtils.subModulo(x, y, N);
				tmpBSum += SFTUtils.innerProduct(SFTUtils.chi(N,v,y),query.get(x_sub_y));
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
}
