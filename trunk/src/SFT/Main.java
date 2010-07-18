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
		// testing on wav file - direct product function over Z_N
		test1();
	}
	
	/*
	 * Testing on wav file: direct product function over Z_N
	 */
	private static void test1() throws Exception{
		long[] G = new long[]{30000};
		long ms = 2*SFTUtils.calcLogG(G)[0];
		int numOfIterations = 2;
		double tau = 0.015;
		
		// create function from WAV file (using MATLAB output from WAV to CSV)
		DirectProdFunction wavFunc = new WavFunction("matlab\\wav\\output.csv",G);
		
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		SFTUtils.ResultFunction f_tag = new SFTUtils.ResultFunction(G,
				SFT.getSignificantElements(G, tau, wavFunc, ms, ms, numOfIterations));
		
		/*
		// writing the result function into a binary file
		FileOutputStream output = new FileOutputStream("matlab\\wav\\java_output_func_2iters.bin");
		ObjectOutputStream out = new ObjectOutputStream(output);
		out.writeObject(f_tag);
		out.close();
		System.out.println(">>> wrote function to binary file");
		*/
		
		// write result function into txt file to be parsed into WAV in MATLAB
		FileWriter f = new FileWriter("matlab\\wav\\java_output_2iters.txt");
		PrintWriter p = new PrintWriter(f);
		for(int i=0; i<G[0]; i++){
			Complex val = f_tag.getValue(new long[]{i});
			p.print(val.getRe()+" "+val.getIm()+"\n");
		}
		p.close(); f.close();
	}
}
