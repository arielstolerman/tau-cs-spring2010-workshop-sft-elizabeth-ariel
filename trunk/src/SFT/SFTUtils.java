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

import java.util.*;
import Function.DirectProdFunction;
import Function.FiniteAbelianFunction;
import Function.FunctionException;

public class SFTUtils {

	
	/* *****************************************
	 * General calculations and parameters check
	 *******************************************/
	
	/**
	 * @param elem	vector to print in format (x1,...,xk)
	 * @return		string representation of the vector
	 */
	public static String vectorToString(long[] elem){
		String ans = "(";
		int k = elem.length;
		for (int i=0; i<k; i++){
			ans += elem[i]+",";
		}
		ans = ans.substring(0, ans.length()-1)+")";
		return ans;
	}
	
	public static String vectorToString(Long[] Elem){
		long[] elem = new long[Elem.length];
		for(int i=0; i<Elem.length; i++) elem[i] = Elem[i].longValue();
		return vectorToString(elem);
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
	protected static double calcDelta(double delta_t, double coeff, double fEucNorm, double tau, long N){
		return delta_t/( coeff * Math.pow(Math.pow(fEucNorm,2)/tau,1.5) *
				(Math.log(N)/Math.log(2)) );
	}
	
	/**
	 * parameters check - default parameters
	 * @param G
	 * @param delta
	 * @param tau
	 * @throws SFTException
	 */
	protected static void checkParameters(long[] G, double delta, double tau, double fInfNorm,
			double fEuclideanNorm, float deltaCoeff, float randSetsCoeff) throws SFTException{
		if (G == null || G.length < 1){
			throw new SFTException("G must be bigger than 0.");
		}
		for (int i=0; i<G.length; i++){
			if (G[i] <= 0){
				throw new SFTException("all Ns must be positive. G["+i+"] is "+G[i]); 
			}
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
	
	/**
	 * parameters check - norms
	 * @param fInfNorm
	 * @param fEuclideanNorm
	 * @throws SFTException
	 */
	protected static void checkNorms(double fInfNorm, double fEuclideanNorm) throws SFTException{
		
	}
	
	/**
	 * parameters check - constants
	 * @param deltaCoeff
	 * @param randSetsCoeff
	 * @throws SFTException
	 */
	protected static void checkConstants(float deltaCoeff, float randSetsCoeff) throws SFTException{
		
	}
	
	/* ***************************
	 * finite Abelian adaptation
	 *****************************/
	
	/**
	 * Calculates L from L' (for the finite Abelian implementation)
	 * @param Ltag
	 * @param G
	 * @return
	 */
	protected static Set<Long> getAbelianRepresentation(Set<long[]> Ltag, long[][] G){
		Set<Long> L = new HashSet<Long>();
		for (long[] elem: Ltag){
			Long newElem = calcAbelianProd(elem, G);
			L.add(newElem);
		}
		return L;
	}
	
	protected static long[] getGFromAbelianFunc(FiniteAbelianFunction func){
		long[][] funcG = func.getG();
		long[] G = new long[funcG.length];
		for (int i=0; i<G.length; i++) G[i] = funcG[i][1];
		return G;
	}
	
	protected static Long[] getDirectProdGFromAbelianG(Long[][] abelianG){
		Long[] G = new Long[abelianG.length];
		for (int i=0; i<G.length; i++) G[i] = abelianG[i][1];
		return G;
	}
	
	/**
	 * Isomorphism from an element in a direct product G to an element in a finite Abelian G 
	 * @param elem
	 * @param G
	 * @return
	 */
	protected static Long calcAbelianProd(long[] elem, long[][] G){
		long newElem = 1;
		for (int i=0; i<elem.length; i++){
			long g = G[i][0];
			long N = G[i][1];
			long x = elem[i];
			newElem *= (g*x % N);
		}
		return new Long(newElem);
	}
	
	protected static Long calcAbelianProd(long[] elem, Long[][] G){
		long[][] smallG = new long[G.length][];
		for (int i=0; i<G.length; i++) smallG[i] = new long[G[i].length];
		for (int i=0; i<G.length; i++)
			for (int j=0; j<G[i].length; j++)
				smallG[i][j] = G[i][j].longValue();
		return calcAbelianProd(elem,smallG);
	}
	
	protected static Long calcAbelianProd(Long[] elem, Long[][] G){
		long[] e = new long[elem.length];
		for (int i=0; i<e.length; i++) e[i] = elem[i].longValue();
		return calcAbelianProd(e,G);
	}
	
	/**
	 * class for disguising a FiniteAbelianFunction as a DirectProdFunction
	 */
	protected static class DirectedProdFromAbelianFunc extends DirectProdFunction{
		protected FiniteAbelianFunction func;
		protected long[][] funcG;
		
		public DirectedProdFromAbelianFunc(long[] G, FiniteAbelianFunction func) throws FunctionException {
			super(G);
			this.func = func;
			this.funcG = func.getG();
		}
		
		public Complex getValue(long[] elem){
			long abelianElem = calcAbelianProd(elem, funcG);
			return func.getValue(abelianElem);
		}
		
	}
	
	/* *********************************
	 * Mathematical functions
	 ***********************************/
	
	/**
	 * @param G		an integer array describing the group Z_N1 X ... X Z_Nk
	 * @return		an integer array of log_2(Ni), rounded up, for i in {1,...,k}
	 */
	protected static int[] calcLogG(long[] G){
		int[] res= new int[G.length];
		for (int i=0; i<G.length; i++){
			res[i]=(int)Math.ceil(Math.log(G[i])/Math.log(2));
		}
		return res;
	}
	
	/**
	 * @param N		describing Z_N
	 * @return		log_2(N), rounded up
	 */
	protected static int calcLogN(long N){
		return (int)Math.ceil(Math.log(N)/Math.log(2));
	}
	
	/**
	 * calculate Chi over G
	 * @param t		the size of the vector to look at for the calculation
	 * @param G		vector of values describing G, i.e. Cartesian multiplication of Z_Ni
	 * @param v		the vector of elements in G defining the chi function
	 * @param y		input vector for the chi function
	 * @return		chi_(v)[y] = chi_(alpha_1,...,alpha_k)[y_1,...,y_k]
	 */
	public static Complex chi(int t, long[] G, long[] v, long[] y){
		Complex ans = new Complex(1,0);
		for(int i=0; i<t; i++){
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
	 * @param G		an integer array describing the group Z_N1 X ... X Z_Nk
	 * @return		a set of elements in G, uniformly randomly selected
	 */
	protected static Set<long[]> generateRandomSubsetA(long m_A, long[] G){
		return generateRandomSubset(m_A,G,G.length+1,0);
	}
	
	/**
	 * @param m_B	potential size of the set
	 * @param G		an integer array describing the group Z_N1 X ... X Z_Nk 
	 * @param l		a value between 1 and log(N)
	 * @param t		a value between 1 and k
	 * @return		a set of size m_B, of uniformly randomly selected vectors 
	 * 				over Z_N1 X ... X Z_Nt-1 X {0,...,2^(l-1)-1} X {0} X ... X {0} 
	 * 				(with >= k-t zero coordinates)
	 */
	protected static Set<long[]> generateRandomSubsetBtl(long m_B, long[] G, int t, int l){
		Set<long[]> res = new HashSet<long[]>();
		
		// compatibility with G= Z_N
		if(t==1){
			long pow = (long)Math.pow(2, l-1);
			// if 2^(l-1) < m_B, no need to randomly choose elements for B, take all 0,...,2^(l-1)-1
			if (pow <= m_B){
				// take all elements in {0,...,2^(l-1)-1} to B_tl
				for (long i=0; i<pow; i++){
					long[] elem = new long[G.length];
					elem[0] = i;
					for (int j=1; j<G.length; j++) elem[j]=0;
					res.add(elem);
				}
				return res;
			}
		}
		return generateRandomSubset(m_B,G,t,l);
	}
	
	/**
	 * @param sizeOfSet		the size of the needed set of elements
	 * @param G				an integer array describing the group Z_N1 X ... X Z_Nk
	 * @param l				a value between 1 and log(Nt)a value between 1 and log(N)
	 * @param t				a value between 1 and k
	 * @return				a set of size sizeOfSet, of uniformly randomly selected vectors 
	 * 						over Z_N1 X ... X Z_Nt-1 X {0,...,2^(l-1)-1} X {0} X ... X {0} 
	 * 						(with >= k-t zero coordinates)
	 * 						note that when t=k+1 and l=0, we get a set of vectors over 
	 * 						Z_N1 X ... X Z_Nk
	 */
	protected static Set<long[]> generateRandomSubset(long sizeOfSet, long[] G, int t, long l){
		Set<long[]> res = new HashSet<long[]>();
		long pow;
		int j;
		Random rand = new Random();
		for(long i=0; i<sizeOfSet; i++){
			boolean doAgain;
			long[] e= new long[G.length];
			do{
				for(j=0; j<t-1; j++)
					e[j] = (long)Math.floor(rand.nextDouble()*G[j]);
				if(l>0){
					pow = (long)Math.pow(2, l-1);
					e[j]= (long)Math.floor(rand.nextDouble()*pow);
				for(int k=j+1; k<G.length; k++)
					e[k] = 0;
				}
				if(contains(res, e))
					doAgain = true;
				else
					doAgain = false;
			} while (doAgain);
			res.add(e);
		}	
		return res;
	}
	
	/* ******************
	 * Matlab adjustments
	 * *****************/
	
	/**
	 * A class for holding temporary data between part 1 and part 2 of the SFT algorithm for direct product G.
	 */
	protected static class MatlabTemporaryRepositoryDirectProd{
		private Set<long[]>[][] sets;
		private Long[][] Q;
		private Map<String,Complex> query;
	
		/**
		 * default constructor
		 */
		public MatlabTemporaryRepositoryDirectProd(Set<long[]>[][] sets, Long[][] Q, Map<String,Complex> query){
			this.sets = sets;
			this.Q = Q;
			this.query = query;
		}
		
		// getters
		public Set<long[]>[][] getSets(){return sets;}
		public Long[][] getQ(){return Q;}
		public Map<String,Complex> getQuery(){return query;}
		// setters
		public void setQuery(Map<String,Complex> query){this.query = query;}
	}
	
	/**
	 * A class for holding temporary data between part 1 and part 2 of the SFT algorithm for finite Abelian G.
	 */
	protected static class MatlabTemporaryRepositoryFiniteAbelian{
		private Set<long[]>[][] sets; // stays the same as for direct product
		private Long[] Q;
		private Long[][] directProdQ;
		private Map<String,Complex> query;
	
		/**
		 * default constructor
		 */
		public MatlabTemporaryRepositoryFiniteAbelian(Set<long[]>[][] sets, Long[] Q, Long[][] directProdQ, Map<String,Complex> query){
			this.sets = sets;
			this.Q = Q;
			this.directProdQ = directProdQ;
			this.query = query;
		}
		
		// getters
		public Set<long[]>[][] getSets(){return sets;}
		public Long[] getQ(){return Q;}
		public Long[][] getDirectProdQ(){return directProdQ;}
		public Map<String,Complex> getQuery(){return query;}
		// setters
		public void setQuery(Map<String,Complex> query){this.query = query;}
	}
	
	protected static MatlabTemporaryRepositoryFiniteAbelian getMatlabFiniteAbelianRep(MatlabTemporaryRepositoryDirectProd rep, Long[][] G){
		// generate isomorphic parameters
		// sets:
		Set<long[]>[][] sets = rep.getSets();
		// Q:
		long[][] repQ = new long[rep.getQ().length][];
		for (int i=0; i<rep.getQ().length; i++){
			Long[] tmp = rep.getQ()[i];
			repQ[i] = new long[tmp.length];
			for (int j=0; j<tmp.length; j++) repQ[i][j] = tmp[j].longValue();
		}
		Long[] Q = new Long[repQ.length];
		for (int i=0; i<repQ.length; i++){
			Q[i] = calcAbelianProd(repQ[i],G);
		}
		// directProdQ:
		Long[][] directProdQ = rep.getQ();
		// query
		Map<String,Complex> query = rep.getQuery();
		
		// create new finite Abelian repository and return it
		return new MatlabTemporaryRepositoryFiniteAbelian(sets, Q, directProdQ, query);
	}
	
	protected static MatlabTemporaryRepositoryDirectProd getMatlabDirectProdRep(MatlabTemporaryRepositoryFiniteAbelian rep, Long[][] G){
		// generate isomorphic parameters
		// sets:
		Set<long[]>[][] sets = rep.getSets();
		// Q:
		Long[][] Q = rep.getDirectProdQ();
		
		// query
		Map<String,Complex> newQuery = new HashMap<String,Complex>();
		Map<String,Complex> query = rep.getQuery();
		for(Long[] elem: Q){
			Long abelianElem = calcAbelianProd(elem, G);
			String strElem = vectorToString(elem);
			newQuery.put(strElem, query.get("("+abelianElem.toString()+")"));
		}
		
		// create new direct product repository and return it
		return new MatlabTemporaryRepositoryDirectProd(sets, Q, newQuery);
	}
}
