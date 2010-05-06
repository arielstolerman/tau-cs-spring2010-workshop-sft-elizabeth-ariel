package Function;

import java.util.HashMap;
import java.util.Map;

import SFT.Complex;
import SFT.FunctionException;
import SFT.SFTUtils;

/**
 * Private class Polynomial is used to describe Fourier polynomials over Z_N by elements and their
 * complex coefficients.
 */
public class FourierPolynomial extends Function{
	private Map<String,Complex> terms; // the string representation of the long-vector
	private String id;
	
	/**
	 * default constructor
	 * @param id
	 */
	public FourierPolynomial(long[] G, String id) throws FunctionException{
		super(G);
		this.id = id;
		terms = new HashMap<String,Complex>();
	}
	
	// getters
	
	/**
	 * @param alpha		the element (vector) for which the coefficient is fetched
	 * @return			the coefficient or null if doesn't exist
	 */
	public Complex getCoeff(long[] alpha){
		String vector = SFTUtils.printVector(alpha);
		return terms.get(vector);
	}
	
	/**
	 * @return			the id of the polynomial
	 */
	public String getId() {
		return id;
	}
	
	/**
	 * get the string representation of the polynomial
	 */
	public String toString(){
		String str = "";
		
		for (String elem: terms.keySet()){
			Complex coeff = terms.get(elem);
			str += "("+coeff.toString()+")*chi_"+elem+"[X] + ";
		}
		str = str.substring(0,str.length()-3)+", X in G";
		
		return str;
	}
	
	// setters
	
	/**
	 * adds the term to the polynomial
	 * if already exist, adds the coefficients
	 * @param alpha		the element
	 * @param re		real part
	 * @param im		imaginary part
	 */
	public void addUpdateTerm(long[] alpha, double re, double im){
		Complex coeff = this.getCoeff(alpha);
		if (coeff == null){
			// create new entry
			coeff = new Complex(re,im);
			terms.put(SFTUtils.printVector(alpha),coeff);
		} else {
			// add to existing coefficient
			coeff.addComplex(re, im);
		}
	}
	
	// calculations
	
	/**
	 * @param x:	input for the polynomial p
	 * @return:		the complex value of p(x) which is SUM_(alpha vector in G) [coeff_alpha * chi_alpha-vector(x)]
	 */
	public Complex getValue(long[] x){
		Complex ans = new Complex(0,0);
		
		for(String alpha: terms.keySet()){
			long[] alphaVector = SFTUtils.getVectorFromString(alpha);
			Complex coeff = terms.get(alpha);
			ans.addComplex(Complex.mulComplex(coeff,SFTUtils.chi(this.G,alphaVector,x)));
		}
		
		return ans;
	}
}