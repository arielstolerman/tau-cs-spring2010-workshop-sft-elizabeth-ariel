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

import javax.management.Query;
import javax.media.format.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;

import com.sun.media.format.WavAudioFormat;
import com.sun.media.multiplexer.audio.WAVMux;

import Function.DirectProdFunction;
import Function.Function;
import Function.FunctionException;
import Function.FiniteAbelianFunction;
import Function.FourierPolynomial;
import Function.XMLFourierPolynomial;
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
		//test1();
		//test2();
		//test3();
		//test4();
		
		long[] G = new long[]{1000, 1000, 1000, 1000};
		int numOfIterations = 1;
		double tau = 50000;
		int[] logs = SFTUtils.calcLogG(G);
		long ms = 10;
		String filename = "web\\sample.xml";
		// maximum Q size for given ms: 3,750 which is 3.75E-7 % of G size

		// create function from XML file
		DirectProdFunction p = new XMLFourierPolynomial(new File(filename), G);
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		// run n times and save ONLY the intersection of L from all iterations
		int n = 1;
		int[] sizes = new int[n];
		SFTUtils.ResultFunction f = new SFTUtils.ResultFunction(G,
				SFT.getSignificantElements(G,tau,p,ms,ms,numOfIterations));
		Map<long[],Complex> map = f.getMapping();
		sizes[0] = map.size();

		for (int i=1; i<n; i++){
			f = new SFTUtils.ResultFunction(G,SFT.getSignificantElements(
					G,tau,p,ms,ms,numOfIterations));
			Map<long[],Complex> tempMap = f.getMapping();

			Set<long[]> toRemove = new HashSet<long[]>();
			for(long[] elem: map.keySet()){
				String e = SFTUtils.vectorToString(elem);
				boolean removeElem = true;
				for (long[] t: tempMap.keySet()){
					if (SFTUtils.vectorToString(t).equals(e)){
						removeElem = false;
						break;
					}
				}
				if (removeElem) toRemove.add(elem);
			}
			for(long[] elem: toRemove) map.remove(elem);
			sizes[i] = map.size();
		}

		System.out.println("Final mapping after "+n+" runs of SFT:");
		for(long[] elem: map.keySet()){
			System.out.println(SFTUtils.vectorToString(elem)+": "+map.get(elem));
		}
		System.out.println("sizes:");
		for(int i=0; i<n; i++){
			System.out.println("sizes["+i+"]: "+sizes[i]);
		}
	}
	
	/*
	 * Testing on wav file: direct product function over Z_N
	 */
	private static void test1() throws Exception{
		long[] G = new long[]{30000};
		long ms = 30;
		int numOfIterations = 1;
		double tau = 0.04;
		String filename = "matlab\\FINAL\\samples\\orchestra_short\\orchestra_short";
		String longFileName = filename+"_tau-"+(tau)+"_ma-"+ms+"_mb-"+ms+"_iters-"+numOfIterations+"_recs-1_FIX";
		
		// create function from WAV file (using MATLAB output from WAV to CSV)
		DirectProdFunction wavFunc = new WavFunction(filename+".csv",G);
		
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		SFTUtils.ResultFunction f_tag = new SFTUtils.ResultFunction(G,
				SFT.getSignificantElements(G, tau, wavFunc, ms, ms, numOfIterations));
		
		// write result function into txt file to be parsed into WAV in MATLAB
		FileWriter f = new FileWriter(longFileName+".txt");
		PrintWriter p = new PrintWriter(f);
		for(int i=0; i<G[0]; i++){
			Complex val = f_tag.getValue(new long[]{i});
			p.print(val.getRe()+" "+val.getIm()+"\n");
		}
		p.close(); f.close();
	}
	
	/*
	 * Testing on wav file using Matlab's methods
	 */
	private static void test2() throws Exception{
		Long[] G = new Long[]{new Long(30000)};
		long[] g = new long[]{30000};
		long ms = SFTUtils.calcLogG(g)[0];
		int numOfIterations = 1;
		double tau = 0.1;
		String filename = "sample06_short";
		String longFileName = filename+"_ms-"+(ms)+"_iters-"+(numOfIterations)+"_tau-"+(tau)+"_matlab";
		
		// create function from WAV file (using MATLAB output from WAV to CSV)
		DirectProdFunction wavFunc = new WavFunction("matlab\\wav\\"+filename+".csv",g);
		
		// RUN THE SFT ALGORITHM (<<< MATLAB VERSION >>>) TO APPROXIMATE THE FUNCTION
		// run part 1
		MatlabTemporaryRepositoryDirectProd rep = SFT.runMatlabSFTPart1Internal(new Boolean(true),G,tau,ms,ms);
		
		// calculate the function
		Map<String,Complex> query = new HashMap<String,Complex>();
		System.out.println(">>>>>>>>>>>>> creating query");
		for(long i=0; i<G[0]; i++) query.put("("+i+")", wavFunc.getValue(new long[]{i}));
		System.out.println(">>>>>>>>>>>>> done creating query");
		rep.setQuery(query);
		
		// run part 2
		SFTUtils.MatlabTemporaryResultDirectProd res = SFT.runMatlabSFTPart2Internal(G, tau, rep, 1);
		
		Map<long[],Complex> mapping = new HashMap<long[],Complex>();
		Long[][] keys = res.getKeys();
		Long[][] randSet = res.getRandSet();
		int size = G.length;
		for(int i=0; i<keys.length; i++){
			long[] elem = new long[size];
			for(int j=0; j<size; j++) elem[j] = keys[i][j].longValue();
			//Complex val = new Complex(vals[i][0],vals[i][1]);
			//mapping.put(elem, val);
		}
		
		SFTUtils.ResultFunction f_tag = new SFTUtils.ResultFunction(g,mapping);
		
		// write result function into txt file to be parsed into WAV in MATLAB
		FileWriter f = new FileWriter("matlab\\wav\\"+longFileName+"_java.txt");
		PrintWriter p = new PrintWriter(f);
		for(int i=0; i<G[0]; i++){
			Complex val = f_tag.getValue(new long[]{i});
			p.print(val.getRe()+" "+val.getIm()+"\n");
		}
		p.close(); f.close();
	}
	
	private static void test3() throws Exception{
		long[] G = new long[]{30000};
		int numOfIterations = 1;
		double tau = 0.01;
		double infNorm = 1.33970418549327;
		double eucNorm = 0.46955366451643;
		double delta_t = 0.2;
		float deltaCoeff = (float)0.001;
		float maCoeff = (float)0.0000002;
		float mbCoeff = (float)0.0000002;
		float etaCoeff = (float)1;
		String filename = "sample06_short";
		String longFileName = filename+"_delta-"+(delta_t)+"_iters-"+(numOfIterations)+"_tau-"+(tau)+
		"_d-coeff-"+deltaCoeff+"_maCoeff-"+maCoeff+"_mbCoeff-"+mbCoeff+"_etaCoeff-"+etaCoeff;
		
		// create function from WAV file (using MATLAB output from WAV to CSV)
		DirectProdFunction wavFunc = new WavFunction("matlab\\wav\\"+filename+".csv",G);
		
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		SFTUtils.ResultFunction f_tag = new SFTUtils.ResultFunction(G,
				SFT.getSignificantElements(G, tau, wavFunc, numOfIterations, delta_t,
						infNorm, eucNorm, deltaCoeff, maCoeff, mbCoeff, etaCoeff));
		
		// write result function into txt file to be parsed into WAV in MATLAB
		FileWriter f = new FileWriter("matlab\\wav\\"+longFileName+"_java.txt");
		PrintWriter p = new PrintWriter(f);
		for(int i=0; i<G[0]; i++){
			Complex val = f_tag.getValue(new long[]{i});
			p.print(val.getRe()+" "+val.getIm()+"\n");
		}
		p.close(); f.close();
	}
	
	// XML function
	public static void test4() throws Exception{
		long[] G = new long[]{1000,1000};
		int numOfIterations = 1;
		double tau = 60000;
		int[] logs = SFTUtils.calcLogG(G);
		long ms = (long)(0.3*logs[0]*logs[1]);
		String filename = "sample";
		
		// create function from XML file
		DirectProdFunction p = new XMLFourierPolynomial(new File("web\\"+filename+".xml"), G);
		
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		// run n times and save ONLY the intersection of L from all iterations
		int n = 10;
		int[] sizes = new int[n];
		SFTUtils.ResultFunction f = new SFTUtils.ResultFunction(G,
				SFT.getSignificantElements(G,tau,p,ms,ms,numOfIterations));
		Map<long[],Complex> map = f.getMapping();
		sizes[0] = map.size(); 
			
		for (int i=1; i<n; i++){
			f = new SFTUtils.ResultFunction(G,SFT.getSignificantElements(G,tau,p,ms,ms,numOfIterations));
			Map<long[],Complex> tempMap = f.getMapping();
			
			Set<long[]> toRemove = new HashSet<long[]>();
			for(long[] elem: map.keySet()){
				String e = SFTUtils.vectorToString(elem);
				boolean removeElem = true;
				for (long[] t: tempMap.keySet()){
					if (SFTUtils.vectorToString(t).equals(e)){
						removeElem = false;
						break;
					}
				}
				if (removeElem) toRemove.add(elem);
			}
			for(long[] elem: toRemove) map.remove(elem);
			sizes[i] = map.size();
		}
		
		System.out.println("Final mapping after "+n+" runs of SFT:");
		for(long[] elem: map.keySet()){
			System.out.println(SFTUtils.vectorToString(elem)+": "+map.get(elem));
		}
		System.out.println("sizes:");
		for(int i=0; i<n; i++){
			System.out.println("sizes["+i+"]: "+sizes[i]);
		}
	}
	
	/**
	 * write function into binary file named "<filename>.bin"
	 * @param func
	 * @param filename
	 */
	private static void funcToBin(SFTUtils.ResultFunction func, String filename){
		try{
		// writing the result function into a binary file
		FileOutputStream output = new FileOutputStream("matlab\\wav\\"+filename+".bin");
		ObjectOutputStream out = new ObjectOutputStream(output);
		out.writeObject(func);
		out.close();
		System.out.println(">>> wrote function to binary file.");
		} catch (IOException ioe){
			System.err.println(">>> could not write function to binary file.");
		}
	}
}
