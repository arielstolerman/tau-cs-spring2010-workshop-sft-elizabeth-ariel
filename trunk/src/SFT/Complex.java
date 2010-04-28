/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: Complex.java
 * Description: Code for class Complex, implementation of complexed numbers
 * 
 */

package SFT;

/**
 * @author Elizabeth Firman and Ariel Stolerman
 * This class is an implementaiton for representing Complex numbers.
 */
public class Complex {
	// class members
	private double real, imaginary;
	
	/*
	 * constructors
	 */
	
	/**
	 * Constructs a new Complex object with real and imaginary coordinates.
	 * @param real
	 * 				The real coordinate
	 * @param imaginary
	 * 				The imaginary coordinate
	 */
	public Complex(double real, double imaginary){
		this.real = real;
		this.imaginary = imaginary;
	}
	
	/*
	 * methods
	 */
	
	// getters
	
	/**
	 * Returns the real coordinate of the complexed number.
	 * @return
	 * 			The real coordinate of the complexed number.
	 */
	public double getRe(){
		return this.real;
	}
	
	/**
	 * Returns the imaginary coordinate of the complexed number.
	 * @return
	 * 			The imaginary coordinate of the complexed number.
	 */
	public double getIm(){
		return this.imaginary;
	}
	
	/**
	 * Returns an array of two doubles, where the first element is the number's real coordinate
	 * and the second element is the number's imaginary coordinate.
	 * @return
	 * 			A double array with the real coordinate as the first element
	 * 			and the imaginary coordinate as the second element.
	 */
	public double[] getComplexCoords(){
		double[] res = new double[2];
		res[0] = this.real;
		res[1] = this.imaginary;
		
		return res;
	}
	
	/**
	 * Returns the string representation of the complex number in the form "a + bi".
	 * @return
	 * 			The string representation of the complex number in the form "a + bi".
	 */
	public String toString(){
		return real+(imaginary<0?"":"+")+imaginary+"i";
	}
	
	/**
	 * Calculates and returns the Euclidean norm of the complexed number
	 * @return
	 * 			The Euclidean norm of the complexed number
	 */
	public double getEuclideanNorm(){
		return Math.pow(real, 2)+Math.pow(imaginary, 2);
	}
	
	// setters
	
	/**
	 * Sets the real and imaginary coordinates of the number to the new given coordinates.
	 * @param real
	 * 				The new real coordinate.
	 * @param imaginary
	 * 				The new imaginary coordinate.
	 */
	public void setComplex(double real, double imaginary){
		this.real = real;
		this.imaginary = imaginary;
	}
	
	/**
	 * Adds the given complex number to this one.
	 * @param real
	 * 				The real value to be added.
	 * @param imaginary
	 * 				The imaginary value to be added.
	 */
	public void addComplex(double real, double imaginary){
		this.real = this.real+real;
		this.imaginary = this.imaginary+imaginary;
	}
	
	/**
	 * Adds the given complex number to this one.
	 * @param complex
	 * 				The complexed number to be added.
	 */
	public void addComplex(Complex complex){
		this.addComplex(complex.getRe(), complex.getIm());
	}
	
	// static mathematical functions
	
	/**
	 * Returns the multiplication of the two given complexed numbers.
	 * @param complex1
	 * 					The first complexed number for the multiplication.
	 * @param complex2
	 * 					The second complexed number for the multiplication.
	 * @return
	 * 					The multiplication of the two given complexed numbers.
	 */
	public static Complex mulComplex(Complex complex1, Complex complex2){
		// (a + bi)(c + di) = (ac - bd) + (bc + ad)i
		return new Complex(
				complex1.getRe()*complex2.getRe() - complex1.getIm()*complex2.getIm(),
				complex1.getIm()*complex2.getRe() + complex1.getRe()*complex2.getIm()
				); 
	}
}
