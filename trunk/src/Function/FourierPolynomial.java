/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: FourierPolynomial.java
 * Description: Code for class FourierPolynomial, implementation of the abstract class Function for Fourier polynomials
 * 
 */

package Function;

import java.util.*;
import SFT.*;

/**
 * Describes a Fourier polynomial
 * by a list of terms and their coefficients, i.e. for a Fourier polynomial:<br>
 * p(x) = &sum;c<sub>&alpha;</sub>&bull;&Chi;<sub>&alpha;</sub>(x) it holds the mapping of &alpha;
 * to its coefficient c<sub>&alpha;</sub>.
 * @author Elizabeth Firman and Ariel Stolerman
 */
public class FourierPolynomial extends DirectProdFunction{
	private Map<String,Complex> terms; // the string representation of the long-vector
	private String id;
	
	/**
	 * Default constructor.
	 * @param id
	 * 			unique identifier for the polynomial.
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
		String vector = SFTUtils.vectorToString(alpha);
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
	
	/**
	 * get the function as a Matlab script
	 */
	public String toMatlabScript(String scriptName){
		String script = "function[res]="+scriptName+"(x,G)\n";
		script += "% this is an automatically generated function from a FourierPolynomial java object.\n\n";
		script += "k = length(G);\n";
		script += "alpha = [";
		// insert all elements, each element in a column
		String[] elemsArr = new String[terms.keySet().size()];
		terms.keySet().toArray(elemsArr);
		for (int i=0; i<elemsArr.length; i++){
			String vec = elemsArr[i].replace(',', ' ');
			vec = vec.substring(1, vec.length()-1);
			script += vec+";";
		}
		script = script.substring(0, script.length()-1)+"];\n";
		script += "alpha = transpose(alpha);\n";
		// insert all coeffs
		script += "coeff_alpha = [";
		for (int i=0; i<elemsArr.length; i++){
			script += "("+terms.get(elemsArr[i]).toString()+") ";
		}
		script += "];\n\n";
		
		script += "res = 0;\n";
		script += "s = size(alpha); alpha_size = s(1)*s(2);\n";
		script += "for j=1:k:alpha_size;\n";
		script += "\t"+"vec = alpha(j:(j+k-1));\n";
		script += "\t"+"tmp = 1;\n";
		script += "\t"+"for l=1:k;\n";
		script += "\t\t"+"curr_alpha = vec(l);\n";
		script += "\t\t"+"term = 2*pi*curr_alpha*x(l)/G(l);\n";
		script += "\t\t"+"tmp = tmp * exp(term*i);\n";
		script += "\t"+"end\n";
		script += "\t"+"coeff_index = floor(j/k)+1;\n";
		script += "\t"+"res = res + coeff_alpha(coeff_index)*tmp;\n";
		script += "end\n";
		
		return script;
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
			terms.put(SFTUtils.vectorToString(alpha),coeff);
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
			ans.addComplex(Complex.mulComplex(coeff,SFTUtils.chi(G.length,this.G,alphaVector,x)));
		}
		
		return ans;
	}
}