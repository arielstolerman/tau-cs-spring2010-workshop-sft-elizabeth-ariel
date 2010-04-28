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

	
	/* *********************************
	 * General calculations
	 ***********************************/
	
	/**
	 * @param	all parameters needed to calculate delta
	 * @return	delta
	 */
	protected static double calcDelta(double delta_t,float coeff, double fEucNorm, double tau, long N){
		return delta_t/( coeff * Math.pow(Math.pow(fEucNorm,2)/tau,1.5) *
				(Math.log(N)/Math.log(2)) );
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
	 * calculate Chi
	 * @param N		describing Z_N
	 * @param v		floor( (a+b)/2 )
	 * @param y		input for the chi function
	 * @return		chi_(floor[(a+b)/2]) (y)
	 */
	protected static Complex chi(long N, long v, long y){
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
	 * subtraction modulo N
	 * @param a
	 * @param b
	 * @param N
	 * @return
	 */
	protected static long subModulo(long a, long b, long N){
		long ans = a - b;
		if (ans < 0) ans += N;
		return ans;
	}
	
	/* *********************************
	 * Random subsets generation
	 ***********************************/
	
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
