package org.coh.mattiyeh.datacruncher.genome;

public class SbsStats {

	private int numCToT = 0;
	private int numCToA = 0;
	private int numCToG = 0;
	private int numTToC = 0;
	private int numTToA = 0;
	private int numTToG = 0;

	public void incrementCToT() {
		numCToT++;
	}

	public void incrementCToA() {
		numCToA++;
	}

	public void incrementCToG() {
		numCToG++;
	}

	public void incrementTToC() {
		numTToC++;
	}

	public void incrementTToA() {
		numTToA++;
	}

	public void incrementTToG() {
		numTToG++;
	}

	public int getNumCToT() {
		return numCToT;
	}

	public int getNumCToA() {
		return numCToA;
	}

	public int getNumCToG() {
		return numCToG;
	}

	public int getNumTToC() {
		return numTToC;
	}

	public int getNumTToA() {
		return numTToA;
	}

	public int getNumTToG() {
		return numTToG;
	}

}
