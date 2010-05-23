package SFT;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map;
import java.util.Set;
import Function.*;
import SFT.*;
import SFT.SFTUtils.MatlabTemporaryRepository;

/*
 * Main for debugging
 */

public class Main {
	/**
	 * main for debugging
	 * @param args
	 */
	public static void main(String[] args) {
		File xmlInput = new File("d:\\sft_tmp\\test2.xml");
		//File xmlInput = new File("d:\\tmp\\test.xml");
		//long[] G = new long[]{Long.valueOf("100000"),Long.valueOf("100000")};
		long[] G = new long[]{Long.valueOf("10000000000")};
		Long[] bigG = new Long[]{Long.valueOf("10000000000")};
		try {
			MatlabTemporaryRepository matlabRep =
				SFT.runMatlabSFTPart1Internal(bigG,0.01,200,(double)28.41,(double)20.0,(float)1,(float)0.0001,new Boolean(true));
			
			Function poly = new XMLFourierPolynomial(xmlInput, G);
			HashMap<String,Complex> query = new HashMap<String,Complex>();
			for(Long[] Elem: matlabRep.getQ()){
				long[] elem = new long[Elem.length];
				for(int i=0; i<Elem.length; i++) elem[i] = Elem[i].longValue();
				query.put(SFTUtils.vectorToString(elem),poly.getValue(elem));
			}
			matlabRep.setQuery(query);
			
			Long[][] L =
				SFT.runMatlabSFTPart2Internal(bigG, 200, matlabRep);
			
		} catch (SFTException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (FunctionException fe){
			Debug.log(">>> FunctionException thrown: "+fe.getMessage());
		}
		/*
		try{
			// get polynomial
			Function poly = new XMLFourierPolynomial(xmlInput, G);
			// calculate the function and output values
			// print the function
			System.out.println("list of polynomials:");
			for (FourierPolynomial p: ((XMLFourierPolynomial)poly).getPolynomials().values()){
				System.out.println(">> "+p.toString());
			}
			//System.out.println("infinity norm: "+poly.calcInfinityNorm());
			//System.out.println("Euclidean norm: "+poly.calcEuclideanNorm());
			
			Set<long[]> L = SFT.runMainSFTAlgorithm(G, 0.01, 200, poly, (double)28.41, (double)20.0, (float)1, (float)0.0001);
			String res = "";
			for(long[] e: L){
				res += SFTUtils.vectorToString(e)+" ";
			}
			System.out.println("\nSFT: "+res+"\n\nDone!");
			
		} catch (FunctionException fe){
			Debug.log(">>> FunctionException thrown: "+fe.getMessage());
		} catch (SFTException se){
			Debug.log(">>> SFTException thrown: "+se.getMessage());
		}
		
		/* **********************************************************
		 * - generate random functions f
		 * - run the algorithm to create the approximated function f'
		 * - check the error percentage of f'-f
		 ************************************************************/
		/*
		// testing constants and variables
		int NUM_OF_POLYS = 1;
		double NON_ZEROS_PERCENTAGE = 0.01;
		double F_VALUES_BOUND = 100; // will bound the function's coeffs between +-F_VALUES_BOUND/2
		
		long N_TEST = 10000;
		double DELTA = 0.1;
		double TAU = 10000;
		float deltaCoeff = 1;
		float randSetsCoeff = (float)0.000001;
		
		FourierPolynomial[] polys = new FourierPolynomial[NUM_OF_POLYS];
		
		// test
		int i=0;
		// create Fourier polynomials
		for (; i<NUM_OF_POLYS; i++){
			try {
				polys[i] = new FourierPolynomial(N_TEST,i+"");
			} catch (FunctionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			// randomly pick elements in Z_N and their coefficients
			long numOfElements = (long)Math.floor(N_TEST*NON_ZEROS_PERCENTAGE);
			System.out.println("num of non zero elements is "+numOfElements);
			for (int j=0; j<numOfElements; j++){
				long alpha = (long)Math.floor(Math.random()*N_TEST);
				polys[i].addUpdateTerm(alpha, getRand(F_VALUES_BOUND), getRand(F_VALUES_BOUND));
			}
			System.out.println("the poly is: "+polys[i].toString());
		}
		
		Debug.log(">>> MAIN TEST: created fourier polynomials >>>");
		
		// calculate significant elements for each polynomial and create it's corresponding polynomial
		FourierPolynomial[] SFTPolys = new FourierPolynomial[NUM_OF_POLYS];
		for(i=0; i<NUM_OF_POLYS;i++){
			try {
				SFTPolys[i] = new FourierPolynomial(N_TEST,i+"");
				Debug.log(">>> MAIN TEST: initialized guess polynomial >>>");
			} catch (FunctionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			double infNorm = polys[i].calcInfinityNorm();
			Debug.log(">>> MAIN TEST: calculated infinity norm >>>");
			double eucNorm = polys[i].calcEuclideanNorm();
			Debug.log(">>> MAIN TEST: calculated infinityEuclidean norm >>>");
			
			// run the SFT algorithm and create polynomial
			try {
				//Debug.setDEBUG_MODE(false);
				Set<Long>[] sets = SFT.runMainSFTAlgorithmDividedPart1(N_TEST, DELTA, TAU, infNorm, eucNorm, deltaCoeff, randSetsCoeff);
				Set<Long> Q = sets[sets.length-1];
				Map<Long,Complex> query = new HashMap<Long,Complex>();
				for (long elem: Q){
					query.put(elem, polys[i].getValue(elem));
				}
				Set<Long> alphas = SFT.runMainSFTAlgorithmDividedPart2(N_TEST, TAU, sets, query);
				
				// create guess polynomial
				FourierPolynomial sftPoly = SFTPolys[i];
				for(long alpha:alphas){
					Complex value = polys[i].getValue(alpha);
					sftPoly.addUpdateTerm(alpha, value.getRe(), value.getIm());
				}
			} catch (SFTException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		// compare results
		for (i=0; i<NUM_OF_POLYS; i++){
			// count errors
			long errCount = 0;
			for (long elem=0; elem<N_TEST; elem++){
				Complex fVal = polys[i].getValue(elem);
				Complex fTagVal = SFTPolys[i].getValue(elem);
				if (fVal.getRe()-fTagVal.getRe() != 0 || fVal.getIm()-fTagVal.getIm() != 0)
					errCount++;
			}
			System.out.println(">>> The error percentage is: "+((double)errCount)/((double)N_TEST));
		}
		*/
	}
	
	private static double getRand(double F_VALUES_BOUND){
		return Math.random()*F_VALUES_BOUND-F_VALUES_BOUND/2;
	}
}
