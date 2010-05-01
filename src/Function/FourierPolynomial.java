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
	private Map<Long,Complex> terms;
	private String id;
	
	/**
	 * default constructor
	 * @param id
	 */
	public FourierPolynomial(long N, String id) throws FunctionException{
		super(N);
		this.id = id;
		terms = new HashMap<Long,Complex>();
	}
	
	// getters
	
	/**
	 * @param alpha		the element for which the coefficient is fetched
	 * @return			the coefficient or null if doesn't exist
	 */
	public Complex getCoeff(long alpha){
		for(long elem: terms.keySet()){
			if (elem == alpha)
				return terms.get(elem); // found it, return the coeff
		}
		return null;
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
		
		for (long elem: terms.keySet()){
			Complex coeff = terms.get(elem);
			str += "("+coeff.toString()+")*chi_("+elem+")[x] + ";
		}
		str = str.substring(0,str.length()-3);
		
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
	public void addUpdateTerm(long alpha, double re, double im){
		Complex coeff = this.getCoeff(alpha);
		if (coeff == null){
			// create new entry
			coeff = new Complex(re,im);
			terms.put(alpha,coeff);
		} else {
			// add to existing coefficient
			coeff.addComplex(re, im);
		}
	}
	
	// calculations
	
	/**
	 * @param x:	input for the polynomial p
	 * @return:		the complex value of p(x) which is SUM_(alpha in Z_N) [coeff_alpha * chi_alpha(x)]
	 */
	public Complex getValue(long x){
		Complex ans = new Complex(0,0);
		
		for(long alpha: terms.keySet()){
			Complex coeff = terms.get(alpha);
			ans.addComplex(Complex.mulComplex(coeff,SFTUtils.chi(this.N,alpha, x)));
		}
		
		return ans;
	}
}