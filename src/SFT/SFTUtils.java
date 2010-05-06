/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: SFTUtils.java
 * Description: Code for class SFTUtils, static methods used in the SFT algorithm implementation
 * 
 */

package SFT;

import java.util.HashSet;
import java.util.Random;
import java.util.Set;

public class SFTUtils {

	
	/* *****************************************
	 * General calculations and parameters check
	 *******************************************/
	
	/**
	 * @param elem	vector to print in format (x1,...,xk)
	 * @return		string representation of the vector
	 */
	public static String printVector(long[] elem){
		String ans = "(";
		int k = elem.length;
		for (int i=0; i<k; i++){
			ans += elem[i]+",";
		}
		ans = ans.substring(0, ans.length()-1)+")";
		return ans;
	}
	
	/**
	 * @param vector	the string representation of a vector, as created by printVector
	 * @return			the long array representing the vector
	 */
	public static long[] getVectorFromString(String vector){
		// remove "(", ")"
		vector = vector.substring(1,vector.length()-1);
		String[] elems = vector.split(",");
		long[] res = new long[elems.length];
		for(int i=0; i<res.length; i++){
			res[i] = Long.parseLong(elems[i]);
		}
		return res;
	}
	
	/**
	 * @param	all parameters needed to calculate delta
	 * @return	delta
	 */
	protected static double calcDelta(double delta_t,float coeff, double fEucNorm, double tau, long N){
		return delta_t/( coeff * Math.pow(Math.pow(fEucNorm,2)/tau,1.5) *
				(Math.log(N)/Math.log(2)) );
	}
	
	/**
	 * parameters check
	 * @param N
	 * @param delta
	 * @param tau
	 * @param fInfNorm
	 * @param fEuclideanNorm
	 * @param deltaCoeff
	 * @param randSetsCoeff
	 * @throws SFTException
	 */
	protected static void checkParameters(long N, double delta, double tau,	double fInfNorm,
			double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		if (N <= 0){
			throw new SFTException("N must be positive."); 
		}
		if (delta <= 0 || delta >= 1){
			throw new SFTException("delta must be in (0,1).");
		}
		if (tau <= 0){
			throw new SFTException("tau must be positive.");
		}
		if (fInfNorm < 0){
			throw new SFTException("The infinity norm of the function must be positive.");
		}
		if (fEuclideanNorm < 0){
			throw new SFTException("The Euclidean norm of the function must be positive.");
		}
		if (deltaCoeff <= 0 || randSetsCoeff <= 0){
			throw new SFTException("The coefficients must be positive.");
		}
	}
	
	/* *********************************
	 * Mathematical functions
	 ***********************************/
	
	/**
	 * @param N		describing Z_N
	 * @return		log_2(N), rounded up
	 */
	protected static int calcLogN(long N){
		return (int)Math.ceil(Math.log(N)/Math.log(2));
	}
	
	/**
	 * calculate Chi over G
	 * @param G		vector of values describing G, i.e. Cartesian multiplication of Z_Ni
	 * @param v		the vector of elements in G defining the chi function
	 * @param y		input vector for the chi function
	 * @return		chi_(v)[y] = chi_(alpha_1,...,alpha_k)[y_1,...,y_k]
	 */
	public static Complex chi(long[] G, long[] v, long[] y){
		Complex ans = new Complex(1,1);
		int k = G.length;
		for(int i=0; i<k; i++){
			ans = Complex.mulComplex(ans, chi(G[i],v[i],y[i]));
		}
		return ans;
	}
	
	/**
	 * calculate Chi over Z_N
	 * @param N		describing Z_N
	 * @param v		the element in Z_N defining the chi function
	 * @param y		input for the chi function
	 * @return		chi_(v)[y]
	 */
	public static Complex chi(long N, long v, long y){
		// chi_v (y) = e^(i2pi * v/N * y) = cos(2pi * v/N * y) + i*sin(2pi * v/N * y)
		double arg = 2 * Math.PI * (((double)v)/((double)N)) * ((double)y);
		double re = Math.cos(arg);
		double im = Math.sin(arg);
		
		return new Complex(re, im);
	}
	
	/**
	 * calculate the inner product over C
	 * @param x	first element in the inner product
	 * @param y	first element in the inner product
	 * @return		the inner product <x,y> = sum[x_i * y_i] (i = 1,2)
	 */
	protected static double innerProduct(Complex x, Complex y){
		return x.getRe()*y.getRe()+x.getIm()*y.getIm();
	}
	
	/**
	 * vector subtraction modulo N
	 * @param a
	 * @param b
	 * @param N
	 * @param k
	 * 			vector length
	 * @return
	 */
	protected static long[] subVectorModulo(long[] a, long[] b, long N, int k){
		long[] ans = new long[k];
		for (int i=0; i<k; i++){
			long tmp = a[i] - b[i];
			if (tmp < 0) tmp += N;
			ans[i] = tmp;
		}
		return ans;
	}
	
	/* *********************************
	 * Random subsets generation
	 ***********************************/
	
	/**
	 * @param Q		a set of vectors
	 * @param elem	a vector
	 * @return		true iff Q contains elem
	 */
	protected static boolean contains(Set<long[]> Q, long[] elem){
		int k = elem.length;
		int i;
		for (long[] e: Q){
			for (i=0; i<k; i++){
				if (elem[i] != e[i])
					break;
			}
			if (i == k) return true;
		}
		return false;
	}
	
	/**
	 * @param m_A	the size of the set
	 * @param N		describing Z_N
	 * @return		a set of elements in Z_N, uniformly randomly selected
	 */
	protected static Set<Long> generateRandomSubsetA(long m_A, long N){
		return generateRandomSubset(m_A,N);
	}
	
	/**
	 * @param m_B	potential size of the set
	 * @param N		describing Z_N
	 * @param l		a value between 1 and log(N)
	 * @return		a set of elements in {0,...,2^(l-1)-1}
	 */
	protected static Set<Long> generateRandomSubsetBl(long m_B, long N, int l){
		Set<Long> res;
		
		// if 2^(l-1) < m_B, no need to randomly choose elements for be, take all 0,...,2^(l-1)-1
		long pow = (long)Math.pow(2, l-1);
		if (pow <= m_B){
			res = new HashSet<Long>();
			// take all elements in {0,...,2^(l-1)-1} to B_l
			for (long i=0; i<pow; i++){
				res.add(i);
			}
		}
		// otherwise, choose randomly m_B elements from 0,...,2^(l-1)-1
		else {
			res = generateRandomSubset(m_B,pow);
		}
		
		return res;
	}
	
	/**
	 * @param sizeOfSet	the size of the needed set of elements
	 * @param randBarrier	the barrier for the range of the randomly selected elements
	 * @return				a set of uniformly randomly selected elements in range (0,1,...,randBarrier-1)
	 * 						of size sizeOfSet 
	 */
	protected static Set<Long> generateRandomSubset(long sizeOfSet, long randBarrier){
		Set<Long> res = new HashSet<Long>();
		// assuming sizeOfSet < randBarrier
		
		Random rand = new Random();
		for(long i=0; i<sizeOfSet; i++){
			boolean doAgain;
			long e;
			do{
				e = (long)Math.floor(rand.nextDouble()*randBarrier);
				for (long elem: res){
					if (elem == e){
						doAgain = true;
						break;
					}
				}
				doAgain = false;
			} while (doAgain);
			res.add(e);
		}
		
		return res;
	}
}
