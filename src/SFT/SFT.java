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
	
}
