import SFT.*;
import Function.*;
import java.io.*;
import java.util.*;

public class Test {

	public static void main(String[] args) throws Exception{
		long[] G = new long[]{1000,1000,1000,1000};
		int numOfIterations = 1;
		double tau = 60000;
		int[] logs = SFTUtils.calcLogG(G);
		long ms = (long)(0.001*logs[0]*logs[1]*logs[2]*logs[3]);
		String filename = "sample.xml";
		// maximum Q size for given ms: 3,750 which is 3.75E-7 % of G size

		// create function from XML file
		DirectProdFunction p = new XMLFourierPolynomial(new File(filename), G);
		// RUN THE SFT ALGORITHM TO APPROXIMATE THE FUNCTION
		// run n times and save ONLY the intersection of L from all iterations
		int n = 5;
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
}
