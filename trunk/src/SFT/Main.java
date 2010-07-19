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
		//test1();
		test2();
	}
	
	/*
	 * Testing on wav file: direct product function over Z_N
	 */
	private static void test1() throws Exception{
		long[] G = new long[]{30000};
		long ms = SFTUtils.calcLogG(G)[0];
		int numOfIterations = 10;
		double tau = 0.1;
		String filename = "sample06_short";
		String longFileName = filename+"_ms-"+(ms)+"_iters-"+(numOfIterations)+"_tau-"+(tau);
		
		// create function from WAV file (using MATLAB output from WAV to CSV)
		DirectProdFunction wavFunc = new WavFunction("matlab\\wav\\"+filename+".csv",G);
		
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		SFTUtils.ResultFunction f_tag = new SFTUtils.ResultFunction(G,
				SFT.getSignificantElements(G, tau, wavFunc, ms, ms, numOfIterations));
		
		// write result function into txt file to be parsed into WAV in MATLAB
		FileWriter f = new FileWriter("matlab\\wav\\"+longFileName+"_java.txt");
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
		Complex[] vals = res.getValues();
		int size = G.length;
		for(int i=0; i<keys.length; i++){
			long[] tmp = new long[size];
			for (int j=0; j<size; j++) tmp[j] = keys[i][j].longValue();
			mapping.put(tmp, vals[i]);
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
