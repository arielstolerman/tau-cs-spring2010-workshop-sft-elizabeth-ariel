import SFT.*;
import Function.*;
import java.io.*;
import java.util.*;

public class Test {

	public static void main(String[] args) throws Exception{
		long[] G = new long[]{1000, 1000, 1000, 1000};
		int numOfIterations = 1;
		double tau = 50000;
		long ma = 10;
		long mb = 10;
		String filename = "web\\sample.xml";

		// create function from XML file
		DirectProdFunction p = new XMLFourierPolynomial(new File(filename), G);
		
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		Map<long[],Complex> res = SFT.getSignificantElements(G,tau,p,ma,mb,numOfIterations);
		
		System.out.println("The result of the SFT is:");
		for(long[] elem: res.keySet()){
			System.out.println("\tElement: "+SFTUtils.vectorToString(elem)+"\tCoefficient: "+res.get(elem));
		}
	}
}
