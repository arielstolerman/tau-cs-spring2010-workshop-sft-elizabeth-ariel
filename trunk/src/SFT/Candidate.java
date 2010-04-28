/*
 * 
 * Workshop Learning & Coding Theory Project - The SFT Algorithm
 * TAU, Spring Semester 2010
 * Elizabeth Firman and Ariel Stolerman
 * 
 * Filename: Function.java
 * Description: Code for class Candidate, used in class SFT
 * 
 */

package SFT;

import java.util.*;

public class Candidate {
	
	private Set<long[]> candidates;
	
	/**
	 * default constructor
	 */
	protected Candidate(){
		this.candidates = new HashSet<long[]>();
	}
	
	/**
	 * constructor with initial interval
	 * @param interval:		the initial interval to be added to the candidates group
	 */
	protected Candidate(long[] interval){
		this.candidates = new HashSet<long[]>();
		this.addInterval(interval);
	}
	
	// getters
	
	/**
	 * @return:		the set of intervals for this candidate
	 */
	protected Set<long[]> getSet(){
		return this.candidates;
	}
	
	/**
	 * @param interval:		an interval {a,b} (a,b in 0,...,N)
	 * @return:				returns true iff the interval is contained in this candidate
	 */
	protected boolean belongsTo(long[] interval){
		for (long[] elem: candidates){
			if (elem[0] == interval[0] && elem[1] == interval[1]) return true;
		}
		return false;
	}
	
	// adders and setters
	/**
	 * @param interval:		interval to be added to the candidates group
	 * Assumption:			the interval is not already contained in the candidates group
	 */
	protected void addInterval(long[] interval){
		candidates.add(interval);
	}

}
