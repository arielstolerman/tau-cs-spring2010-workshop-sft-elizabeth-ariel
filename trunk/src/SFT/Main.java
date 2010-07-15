package SFT;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import javax.media.format.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;

import com.sun.media.format.WavAudioFormat;
import com.sun.media.multiplexer.audio.WAVMux;

import Function.*;
import SFT.*;
import SFT.SFTUtils.DirectedProdFromAbelianFunc;
import SFT.SFTUtils.MatlabTemporaryRepositoryDirectProd;
import SFT.SFTUtils.MatlabTemporaryRepositoryFiniteAbelian;
import SFT.SFTUtils.WavFunction;

/*
 * Main for debugging
 */

public class Main {
	/**
	 * main for debugging
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
		/*File wav = new File("matlab\\wav\\sample.wav");
		AudioInputStream ais = AudioSystem.getAudioInputStream(wav);
		System.out.println(ais.getFormat());
		byte[] b = new byte[128];
		String str = "";
		ais.read(b);
		for (int i=0; i<128; i++) str+= b[i]+" ";
		System.out.println(str);
		
		/* ***************
		 * input variables
		 * ***************/
		double delta_t = 0.01;
		/* xml */ /*
		File xmlInput = new File("matlab\\test2.xml");
		long[] G = new long[]{Long.valueOf("10000000000")};
		double tau = 190;
		double infNorm = 28.41; // 286.2467568832943; - for test3.xml
		double eucNorm = 20; // 0.5349658666523577;
		float deltaCoeff = (float)1;
		float randSetsCoeff = (float)0.0001;
		
		/* wav */
		long[] G = new long[]{30000};
		int numOfIterations = 2;
		double tau = 0.015;
		double infNorm = 1.41421356237310; // infNorm = 1.09773977520164;
		double eucNorm = 0.36997028087942; // eucNorm = 0.32740776983111
		
		float deltaCoeff = (float)0.001;
		float randSetsCoeff = (float)0.0000001;
		
		/* *************
		 *  single test
		 * *************/
		
		try {
			DirectProdFunction poly = new WavFunction("matlab\\wav\\output.csv",G);
			
			SFTUtils.ResultFunction f_tag = new SFTUtils.ResultFunction(G,
					SFT.getSignificantElements(G, delta_t, tau, poly, infNorm, eucNorm, deltaCoeff, randSetsCoeff,numOfIterations));
			/*
			FileOutputStream output = new FileOutputStream("matlab\\wav\\java_output_func_2iters.bin");
			ObjectOutputStream out = new ObjectOutputStream(output);
			out.writeObject(f_tag);
			out.close();
			System.out.println(">>> wrote function to binary file");
			*/

			FileWriter f = new FileWriter("matlab\\wav\\java_output_2iters.txt");
			PrintWriter p = new PrintWriter(f);
			for(int i=0; i<G[0]; i++){
				Complex val = f_tag.getValue(new long[]{i});
				p.print(val.getRe()+" "+val.getIm()+"\n");
			}
			p.close(); f.close();
			/*
			SFTUtils.DiffFunction diff = new SFTUtils.DiffFunction(G, poly, f_tag);
			FileWriter f2 = new FileWriter("matlab\\wav\\diff.txt");
			PrintWriter p2 = new PrintWriter(f2);
			for(int i=0; i<WavFunction.WAV_FUNC_FIXED_SIZE; i++){
				Complex val = diff.getValue(new long[]{i});
				p2.print(val.getRe()+" "+val.getIm()+"\n");
			}
			p2.close(); f2.close();
			*/
		} catch (Exception e){
			System.out.println(">>> exception thrown >>>");
			e.printStackTrace();
		}
		
		/* *****************************
		 *  matlab test - direct product
		 * *****************************/ /*
		
		try{
			
			// Matlab DirectProd PART #1
			Long[] bigG = new Long[G.length];
			DirectProdFunction poly = new XMLFourierPolynomial(xmlInput, G);
			for (int i=0; i<G.length; i++) bigG[i] = new Long(G[i]);
			
			// wav sample parameters:
			bigG = new Long[]{new Long(105840)};
			infNorm = 1.25553298056245;
			eucNorm = 0.21671405578766;
			
			MatlabTemporaryRepositoryDirectProd rep = 
				SFT.runMatlabSFTPart1Internal(bigG, delta_t, tau, infNorm, eucNorm, deltaCoeff, randSetsCoeff, new Boolean(true));
			
			// calculate Q's values
			Map<String,Complex> query = new HashMap<String,Complex>();
			for (Long[] elem: rep.getQ()){
				String key = SFTUtils.vectorToString(elem);
				long[] e = SFTUtils.getVectorFromString(key);
				Complex value = poly.getValue(e);
				query.put(key, value);
			}
			rep.setQuery(query);
			
			// Matlab DirectProd PART #2
			Long[][] L = SFT.runMatlabSFTPart2Internal(bigG, tau, rep);
			
		} catch (FunctionException fe){
			Log.log(">>> FunctionException thrown: "+fe.getMessage());
		} catch (SFTException se){
			Log.log(">>> SFTException thrown: "+se.getMessage());
		}
		
		/* *****************************
		 *  matlab test - finite Abelian
		 * *****************************/ /*
		
		try{
			
			// Matlab DirectProd PART #1
			Long[] bigG = new Long[G.length];
			DirectProdFunction poly = new XMLFourierPolynomial(xmlInput, G);
			Long[][] abelianG = new Long[1][2]; abelianG[0][0] = new Long(1); abelianG[0][1] = new Long(G[0]);
			
			for (int i=0; i<G.length; i++) bigG[i] = new Long(G[i]);
			MatlabTemporaryRepositoryFiniteAbelian rep = 
				SFT.runMatlabSFTPart1Internal(abelianG, delta_t, tau, infNorm, eucNorm, deltaCoeff, randSetsCoeff, new Boolean(true));
			Log.log("creating query...");
			// calculate Q's values
			Map<String,Complex> query = new HashMap<String,Complex>();
			for (Long elem: rep.getQ()){
				long[] e = new long[]{elem.longValue()}; 
				String key = SFTUtils.vectorToString(e);
				Complex value = poly.getValue(e);
				query.put(elem.toString(), value);
			}
			rep.setQuery(query);
			Log.log("moving to part 2...");
			
			// Matlab DirectProd PART #2
			Long[] L = SFT.runMatlabSFTPart2Internal(abelianG, tau, rep);
			
		} catch (FunctionException fe){
			Log.log(">>> FunctionException thrown: "+fe.getMessage());
		} catch (SFTException se){
			Log.log(">>> SFTException thrown: "+se.getMessage());
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

		long[] G_TEST = new long[]{10000};
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
				polys[i] = new FourierPolynomial(G_TEST,i+"");
			} catch (FunctionException e) {
				e.printStackTrace();
			}
			// randomly pick elements in Z_N and their coefficients
			long numOfElements = (long)Math.floor(G_TEST*NON_ZEROS_PERCENTAGE);
			System.out.println("num of non zero elements is "+numOfElements);
			for (int j=0; j<numOfElements; j++){
				long alpha = (long)Math.floor(Math.random()*G_TEST);
				polys[i].addUpdateTerm(alpha, getRand(F_VALUES_BOUND), getRand(F_VALUES_BOUND));
			}
			System.out.println("the poly is: "+polys[i].toString());
		}

		Log.log(">>> MAIN TEST: created fourier polynomials >>>");

		// calculate significant elements for each polynomial and create it's corresponding polynomial
		FourierPolynomial[] SFTPolys = new FourierPolynomial[NUM_OF_POLYS];
		for(i=0; i<NUM_OF_POLYS;i++){
			try {
				SFTPolys[i] = new FourierPolynomial(G_TEST,i+"");
				Log.log(">>> MAIN TEST: initialized guess polynomial >>>");
			} catch (FunctionException e) {
				e.printStackTrace();
			}

			double infNorm = polys[i].calcInfinityNorm();
			Log.log(">>> MAIN TEST: calculated infinity norm >>>");
			double eucNorm = polys[i].calcEuclideanNorm();
			Log.log(">>> MAIN TEST: calculated infinityEuclidean norm >>>");

			// run the SFT algorithm and create polynomial
			try {
				//Log.setLogMode(false);
				Set<Long>[] sets = SFT.runMainSFTAlgorithmDividedPart1(G_TEST, DELTA, TAU, infNorm, eucNorm, deltaCoeff, randSetsCoeff);
				Set<Long> Q = sets[sets.length-1];
				Map<Long,Complex> query = new HashMap<Long,Complex>();
				for (long elem: Q){
					query.put(elem, polys[i].getValue(elem));
				}
				Set<Long> alphas = SFT.runMainSFTAlgorithmDividedPart2(G_TEST, TAU, sets, query);

				// create guess polynomial
				FourierPolynomial sftPoly = SFTPolys[i];
				for(long alpha:alphas){
					Complex value = polys[i].getValue(alpha);
					sftPoly.addUpdateTerm(alpha, value.getRe(), value.getIm());
				}
			} catch (SFTException e) {
				e.printStackTrace();
			}
		}

		// compare results
		for (i=0; i<NUM_OF_POLYS; i++){
			// count errors
			long errCount = 0;
			for (long elem=0; elem<G_TEST; elem++){
				Complex fVal = polys[i].getValue(elem);
				Complex fTagVal = SFTPolys[i].getValue(elem);
				if (fVal.getRe()-fTagVal.getRe() != 0 || fVal.getIm()-fTagVal.getIm() != 0)
					errCount++;
			}
			System.out.println(">>> The error percentage is: "+((double)errCount)/((double)G_TEST));
		} */

		/* ********************************************
		 * Testing for finding correct coefficients for
		 * the SFT algorithm
		 **********************************************/

		// method:
		// - create functions with some percentage of the elements as significant (significant coefficients)
		// - run the SFT algorithm on the functions
		// - create sft-guess functions
		// - check that the difference between the original and the sft-guess is ~ delta
		/*
		try {
			// testing variables
			int NUM_OF_POLYS = 1;									// number of functions to generate and test
			long[] G = new long[]{100,100};						// domain specification
			long[] REGULAR_ELEMS_RANGE = new long[]{0,10};			// regular elements coeff range
			long[] SIGNIFICANT_ELEMS_RANGE = new long[]{50,80};		// significant elements rcoeff ange
			double SIGNIFICAT_ELEMS_PERCENTAGE = 0.01;				// the percentage of the elements to be significant
			double NON_ZERO_ELEMS_PERCENTAGE = 0.05;					// number of non-zero coeff elements
			
			// important stuff
			double DELTA_T = 0.1;									// confidence parameter
			double TAU = 200;										// threshold
			float DELTA_COEFF = 1;									// coeff for calculating delta
			float RAND_SETS_COEFF = (float)0.001;							// coeff for creating random subsets 
			
			FourierPolynomial[] polys;
			FourierPolynomial[] sftPolys;

			// initialize variables
			long numOfElements = 1;
			for (int i=0; i<G.length; i++) numOfElements *= G[i];
			
			// STEP 1: create functions randomly
			// ---------------------------------
			polys = new FourierPolynomial[NUM_OF_POLYS];
			for (int i=0; i<NUM_OF_POLYS; i++) polys[i] = new FourierPolynomial(G, i+"");
			
			long numOfNonZeroRegElems = (long)Math.floor(numOfElements*NON_ZERO_ELEMS_PERCENTAGE*(1-SIGNIFICAT_ELEMS_PERCENTAGE));
			System.out.println("number of non-zero regular elemets: "+numOfNonZeroRegElems);
			long numOfNonZeroSigElems = (long)Math.floor(numOfElements*NON_ZERO_ELEMS_PERCENTAGE*SIGNIFICAT_ELEMS_PERCENTAGE);
			System.out.println("number of non-zero SIGNIFICANT elemets: "+numOfNonZeroSigElems);
			
			System.out.println("Polynomials:");
			System.out.println("============");
			// for each function, choose its coefficients randomly
			for (FourierPolynomial poly: polys){				
				// create regular elements
				for (int j=0; j<numOfNonZeroRegElems; j++){
					// create element
					long[] elem = new long[G.length];
					for (int k=0; k<G.length; k++){
						elem[k] = getRandElem(G[k]);
					}
					// insert as new term
					poly.addUpdateTerm(elem, getRandCoeff(REGULAR_ELEMS_RANGE), getRandCoeff(REGULAR_ELEMS_RANGE));
				}
				
				// create significant elements
				// create regular elements
				for (int j=0; j<numOfNonZeroSigElems; j++){
					// create element
					long[] elem = new long[G.length];
					for (int k=0; k<G.length; k++){
						elem[k] = getRandElem(G[k]);
					}
					// insert as new SIGNIFICANT term
					poly.addUpdateTerm(elem, getRandCoeff(SIGNIFICANT_ELEMS_RANGE), getRandCoeff(SIGNIFICANT_ELEMS_RANGE));
				}
				
				// print to screen
				System.out.println("f_"+poly.getId()+"[X] = "+poly);
				System.out.println("\tinfinity norm: "+poly.calcInfinityNorm());
				System.out.println("\tEuclidean norm: "+poly.calcEuclideanNorm());
			}
			
			// STEP 2: create SFT-guess functions
			// ----------------------------------
			sftPolys = new FourierPolynomial[NUM_OF_POLYS];
			for (int i=0; i<NUM_OF_POLYS; i++) sftPolys[i] = new FourierPolynomial(G, i+"_sft");
			
			for (int i=0; i<polys.length; i++){
				FourierPolynomial poly = polys[i];
				double fInfNorm = poly.calcInfinityNorm(); //TODO
				double fEuclideanNorm = poly.calcEuclideanNorm(); //TODO
				Set<long[]> L = SFT.runMainSFTAlgorithm(G, DELTA_T, TAU, poly, fInfNorm, fEuclideanNorm, DELTA_COEFF, RAND_SETS_COEFF);
				FourierPolynomial sftPoly = sftPolys[i];
				
				// fill sftPoly according to results
				for (long[] elem: L){
					Complex coeff = poly.getValue(elem);
					sftPoly.addUpdateTerm(elem, coeff.getRe(), coeff.getIm());
				}
			}
			
			// STEP 3: calculate the difference between f its sft-guess
			// --------------------------------------------------------
			// create matlab functions from the functions and the sft-functions
			for (int i=0; i<polys.length; i++){
				String f_name = "matlab\\func_"+polys[i].getId()+".m";
				String f_sft_name = "matlab\\func_"+sftPolys[i].getId()+".m";
				
				BufferedWriter f = new BufferedWriter(new FileWriter(f_name));
				BufferedWriter f_sft = new BufferedWriter(new FileWriter(f_sft_name));
				
				f.write(polys[i].toMatlabScript(f_name));
				f_sft.write(sftPolys[i].toMatlabScript(f_sft_name));
			}
			
		} catch (FunctionException e) {
			System.err.println("FunctionException thrown");
		} catch (SFTException e) {
			System.err.println("SFTException thrown");
		} catch (IOException e) {
			e.printStackTrace();
		}
		/* */
	}
	
	private static double getRandCoeff(long[] range){
		return (Math.random()*(range[1]-range[0])+range[0]);
	}
	
	private static long getRandElem(long range){
		return (long)Math.floor(Math.random()*range);
	}

	/*
	private static double getRand(double F_VALUES_BOUND){
		return Math.random()*F_VALUES_BOUND-F_VALUES_BOUND/2;
	} */
}
