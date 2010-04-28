package SFT;

import java.io.File;
import java.util.Set;

import Function.Function;
import Function.XMLFourierPolynomial;

/*
 * Main for debugging
 */

public class Main {
	
	/**
	 * main for debugging
	 * @param args
	 */
	public static void main(String[] args) {
		File xmlInput = new File("d:\\tmp\\test.xml");
		long N = Long.valueOf("10000000000"); 
		
		try{
			// get polynomial
			Function poly = new XMLFourierPolynomial(xmlInput, N);
			// calculate the function and output values
			Set<Long> res = SFT.getSignificatElements(N, 0.1, 200.0, poly, 28.41, 20.0, (float)1.0, (float)0.0001);
			System.out.println("The significat elements are:");
			for (long e: res){
				System.out.println(e+" ");
			}
			
		} catch (FunctionException fe){
			Debug.log(">>> FunctionException thrown: "+fe.getMessage());
		} catch (SFTException se){
			Debug.log(">>> SFTException thrown: "+se.getMessage());
		}
	}
}
