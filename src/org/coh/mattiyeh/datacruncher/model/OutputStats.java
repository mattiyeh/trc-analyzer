package org.coh.mattiyeh.datacruncher.model;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

public class OutputStats {

	private int numSpecimens;
	private int numSamples;
	private int numSpecimensWithMutationData;
	private int numSpecimensWithExpressionData;
	private int numSpecimensWithBoth;

	private int numMutations;

	private int numNonPromoterMutations;
	private int numPromoterMutations;

	private int numNonPromoterSbsMutations;
	private int numPromoterSbsMutations;
	private int numPromoterSbsMutationsHighExp;
	private int numPromoterSbsMutationsMidExp;
	private int numPromoterSbsMutationsLowExp;
	private int numPromoterSbsMutationsZeroExp;

	private int numNonPromoterIndelMutations;
	private int numPromoterIndelMutations;
	private int numPromoterIndelMutationsHighExp;
	private int numPromoterIndelMutationsMidExp;
	private int numPromoterIndelMutationsLowExp;
	private int numPromoterIndelMutationsZeroExp;

	private int numNonPromoterMbsMutations;
	private int numPromoterMbsMutations;
	private int numPromoterMbsMutationsHighExp;
	private int numPromoterMbsMutationsMidExp;
	private int numPromoterMbsMutationsLowExp;
	private int numPromoterMbsMutationsZeroExp;

	private int numNonCfsMutations;
	private int numCfsMutations;

	private int numNonCfsSbsMutations;
	private int numCfsSbsMutations;

	private int numNonCfsIndelMutations;
	private int numCfsIndelMutations;

	private int numNonCfsMbsMutations;
	private int numCfsMbsMutations;
	
	public OutputStats() {
		numSpecimens = 0;
		numSamples = 0;
		numSpecimensWithMutationData = 0;
		numSpecimensWithExpressionData = 0;
		numSpecimensWithBoth = 0;

		numMutations = 0;

		numNonPromoterMutations = 0;
		numPromoterMutations = 0;

		numNonPromoterSbsMutations = 0;
		numPromoterSbsMutations = 0;
		numPromoterSbsMutationsHighExp = 0;
		numPromoterSbsMutationsMidExp = 0;
		numPromoterSbsMutationsLowExp = 0;
		numPromoterSbsMutationsZeroExp = 0;

		numNonPromoterIndelMutations = 0;
		numPromoterIndelMutations = 0;
		numPromoterIndelMutationsHighExp = 0;
		numPromoterIndelMutationsMidExp = 0;
		numPromoterIndelMutationsLowExp = 0;
		numPromoterIndelMutationsZeroExp = 0;

		numNonPromoterMbsMutations = 0;
		numPromoterMbsMutations = 0;
		numPromoterMbsMutationsHighExp = 0;
		numPromoterMbsMutationsMidExp = 0;
		numPromoterMbsMutationsLowExp = 0;
		numPromoterMbsMutationsZeroExp = 0;

		numNonCfsMutations = 0;
		numCfsMutations = 0;

		numNonCfsSbsMutations = 0;
		numCfsSbsMutations = 0;

		numNonCfsIndelMutations = 0;
		numCfsIndelMutations = 0;

		numNonCfsMbsMutations = 0;
		numCfsMbsMutations = 0;
	}
	
	public void addNumSpecimens(int numSpecimens) {
		this.numSpecimens += numSpecimens;
	}
	
	public void addNumSamples(int numSamples) {
		this.numSamples += numSamples;
	}
	
	public void addNumSpecimensWithMutationData(int numSpecimensWithMutationData) {
		this.numSpecimensWithMutationData += numSpecimensWithMutationData;
	}
	
	public void addNumSpecimensWithExpressionData(int numSpecimensWithExpressionData) {
		this.numSpecimensWithExpressionData += numSpecimensWithExpressionData;
	}
	
	public void addNumSpecimensWithBoth(int numSpecimensWithBoth) {
		this.numSpecimensWithBoth += numSpecimensWithBoth;
	}
	
	public void addNumMutations(int numMutations) {
		this.numMutations += numMutations;
	}
	
	public void addNumNonPromoterMutations(int numNonPromoterMutations) {
		this.numNonPromoterMutations += numNonPromoterMutations;
	}
	
	public void addNumPromoterMutations(int numPromoterMutations) {
		this.numPromoterMutations += numPromoterMutations;
	}
	
	public void addNumNonPromoterSbsMutations(int numNonPromoterSbsMutations) {
		this.numNonPromoterSbsMutations += numNonPromoterSbsMutations;
	}
	
	public void addNumPromoterSbsMutations(int numPromoterSbsMutations) {
		this.numPromoterSbsMutations += numPromoterSbsMutations;
	}
	
	public void addNumPromoterSbsMutationsHighExp(int numPromoterSbsMutationsHighExp) {
		this.numPromoterSbsMutationsHighExp += numPromoterSbsMutationsHighExp;
	}
	
	public void addNumPromoterSbsMutationsMidExp(int numPromoterSbsMutationsMidExp) {
		this.numPromoterSbsMutationsMidExp += numPromoterSbsMutationsMidExp;
	}
	
	public void addNumPromoterSbsMutationsLowExp(int numPromoterSbsMutationsLowExp) {
		this.numPromoterSbsMutationsLowExp += numPromoterSbsMutationsLowExp;
	}
	
	public void addNumPromoterSbsMutationsZeroExp(int numPromoterSbsMutationsZeroExp) {
		this.numPromoterSbsMutationsZeroExp += numPromoterSbsMutationsZeroExp;
	}
	
	public void addNumNonPromoterIndelMutations(int numNonPromoterIndelMutations) {
		this.numNonPromoterIndelMutations += numNonPromoterIndelMutations;
	}
	
	public void addNumPromoterIndelMutations(int numPromoterIndelMutations) {
		this.numPromoterIndelMutations += numPromoterIndelMutations;
	}
	
	public void addNumPromoterIndelMutationsHighExp(int numPromoterIndelMutationsHighExp) {
		this.numPromoterIndelMutationsHighExp += numPromoterIndelMutationsHighExp;
	}
	
	public void addNumPromoterIndelMutationsMidExp(int numPromoterIndelMutationsMidExp) {
		this.numPromoterIndelMutationsMidExp += numPromoterIndelMutationsMidExp;
	}
	
	public void addNumPromoterIndelMutationsLowExp(int numPromoterIndelMutationsLowExp) {
		this.numPromoterIndelMutationsLowExp += numPromoterIndelMutationsLowExp;
	}
	
	public void addNumPromoterIndelMutationsZeroExp(int numPromoterIndelMutationsZeroExp) {
		this.numPromoterIndelMutationsZeroExp += numPromoterIndelMutationsZeroExp;
	}
	
	public void addNumNonPromoterMbsMutations(int numNonPromoterMbsMutations) {
		this.numNonPromoterMbsMutations += numNonPromoterMbsMutations;
	}
	
	public void addNumPromoterMbsMutations(int numPromoterMbsMutations) {
		this.numPromoterMbsMutations += numPromoterMbsMutations;
	}
	
	public void addNumPromoterMbsMutationsHighExp(int numPromoterMbsMutationsHighExp) {
		this.numPromoterMbsMutationsHighExp += numPromoterMbsMutationsHighExp;
	}
	
	public void addNumPromoterMbsMutationsMidExp(int numPromoterMbsMutationsMidExp) {
		this.numPromoterMbsMutationsMidExp += numPromoterMbsMutationsMidExp;
	}
	
	public void addNumPromoterMbsMutationsLowExp(int numPromoterMbsMutationsLowExp) {
		this.numPromoterMbsMutationsLowExp += numPromoterMbsMutationsLowExp;
	}
	
	public void addNumPromoterMbsMutationsZeroExp(int numPromoterMbsMutationsZeroExp) {
		this.numPromoterMbsMutationsZeroExp += numPromoterMbsMutationsZeroExp;
	}
	
	public void addNumNonCfsMutations(int numNonCfsMutations) {
		this.numNonCfsMutations += numNonCfsMutations;
	}
	
	public void addNumCfsMutations(int numCfsMutations) {
		this.numCfsMutations += numCfsMutations;
	}
	
	public void addNumNonCfsSbsMutations(int numNonCfsSbsMutations) {
		this.numNonCfsSbsMutations += numNonCfsSbsMutations;
	}
	
	public void addNumCfsSbsMutations(int numCfsSbsMutations) {
		this.numCfsSbsMutations += numCfsSbsMutations;
	}
	
	public void addNumNonCfsIndelMutations(int numNonCfsIndelMutations) {
		this.numNonCfsIndelMutations += numNonCfsIndelMutations;
	}
	
	public void addNumCfsIndelMutations(int numCfsIndelMutations) {
		this.numCfsIndelMutations += numCfsIndelMutations;
	}
	
	public void addNumNonCfsMbsMutations(int numNonCfsMbsMutations) {
		this.numNonCfsMbsMutations += numNonCfsMbsMutations;
	}
	
	public void addNumCfsMbsMutations(int numCfsMbsMutations) {
		this.numCfsMbsMutations += numCfsMbsMutations;
	}

	public int getNumSpecimens() {
		return numSpecimens;
	}

	public void setNumSpecimens(int numSpecimens) {
		this.numSpecimens = numSpecimens;
	}

	public int getNumSamples() {
		return numSamples;
	}

	public void setNumSamples(int numSamples) {
		this.numSamples = numSamples;
	}

	public int getNumSpecimensWithMutationData() {
		return numSpecimensWithMutationData;
	}

	public void setNumSpecimensWithMutationData(int numSpecimensWithMutationData) {
		this.numSpecimensWithMutationData = numSpecimensWithMutationData;
	}

	public int getNumSpecimensWithExpressionData() {
		return numSpecimensWithExpressionData;
	}

	public void setNumSpecimensWithExpressionData(int numSpecimensWithExpressionData) {
		this.numSpecimensWithExpressionData = numSpecimensWithExpressionData;
	}

	public int getNumSpecimensWithBoth() {
		return numSpecimensWithBoth;
	}

	public void setNumSpecimensWithBoth(int numSpecimensWithBoth) {
		this.numSpecimensWithBoth = numSpecimensWithBoth;
	}

	public int getNumMutations() {
		return numMutations;
	}

	public void setNumMutations(int numMutations) {
		this.numMutations = numMutations;
	}

	public int getNumNonPromoterMutations() {
		return numNonPromoterMutations;
	}

	public void setNumNonPromoterMutations(int numNonPromoterMutations) {
		this.numNonPromoterMutations = numNonPromoterMutations;
	}

	public int getNumPromoterMutations() {
		return numPromoterMutations;
	}

	public void setNumPromoterMutations(int numPromoterMutations) {
		this.numPromoterMutations = numPromoterMutations;
	}

	public int getNumNonPromoterSbsMutations() {
		return numNonPromoterSbsMutations;
	}

	public void setNumNonPromoterSbsMutations(int numNonPromoterSbsMutations) {
		this.numNonPromoterSbsMutations = numNonPromoterSbsMutations;
	}

	public int getNumPromoterSbsMutations() {
		return numPromoterSbsMutations;
	}

	public void setNumPromoterSbsMutations(int numPromoterSbsMutations) {
		this.numPromoterSbsMutations = numPromoterSbsMutations;
	}

	public int getNumPromoterSbsMutationsHighExp() {
		return numPromoterSbsMutationsHighExp;
	}

	public void setNumPromoterSbsMutationsHighExp(int numPromoterSbsMutationsHighExp) {
		this.numPromoterSbsMutationsHighExp = numPromoterSbsMutationsHighExp;
	}

	public int getNumPromoterSbsMutationsMidExp() {
		return numPromoterSbsMutationsMidExp;
	}

	public void setNumPromoterSbsMutationsMidExp(int numPromoterSbsMutationsMidExp) {
		this.numPromoterSbsMutationsMidExp = numPromoterSbsMutationsMidExp;
	}

	public int getNumPromoterSbsMutationsLowExp() {
		return numPromoterSbsMutationsLowExp;
	}

	public void setNumPromoterSbsMutationsLowExp(int numPromoterSbsMutationsLowExp) {
		this.numPromoterSbsMutationsLowExp = numPromoterSbsMutationsLowExp;
	}

	public int getNumPromoterSbsMutationsZeroExp() {
		return numPromoterSbsMutationsZeroExp;
	}

	public void setNumPromoterSbsMutationsZeroExp(int numPromoterSbsMutationsZeroExp) {
		this.numPromoterSbsMutationsZeroExp = numPromoterSbsMutationsZeroExp;
	}

	public int getNumNonPromoterIndelMutations() {
		return numNonPromoterIndelMutations;
	}

	public void setNumNonPromoterIndelMutations(int numNonPromoterIndelMutations) {
		this.numNonPromoterIndelMutations = numNonPromoterIndelMutations;
	}

	public int getNumPromoterIndelMutations() {
		return numPromoterIndelMutations;
	}

	public void setNumPromoterIndelMutations(int numPromoterIndelMutations) {
		this.numPromoterIndelMutations = numPromoterIndelMutations;
	}

	public int getNumPromoterIndelMutationsHighExp() {
		return numPromoterIndelMutationsHighExp;
	}

	public void setNumPromoterIndelMutationsHighExp(int numPromoterIndelMutationsHighExp) {
		this.numPromoterIndelMutationsHighExp = numPromoterIndelMutationsHighExp;
	}

	public int getNumPromoterIndelMutationsMidExp() {
		return numPromoterIndelMutationsMidExp;
	}

	public void setNumPromoterIndelMutationsMidExp(int numPromoterIndelMutationsMidExp) {
		this.numPromoterIndelMutationsMidExp = numPromoterIndelMutationsMidExp;
	}

	public int getNumPromoterIndelMutationsLowExp() {
		return numPromoterIndelMutationsLowExp;
	}

	public void setNumPromoterIndelMutationsLowExp(int numPromoterIndelMutationsLowExp) {
		this.numPromoterIndelMutationsLowExp = numPromoterIndelMutationsLowExp;
	}

	public int getNumPromoterIndelMutationsZeroExp() {
		return numPromoterIndelMutationsZeroExp;
	}

	public void setNumPromoterIndelMutationsZeroExp(int numPromoterIndelMutationsZeroExp) {
		this.numPromoterIndelMutationsZeroExp = numPromoterIndelMutationsZeroExp;
	}

	public int getNumNonPromoterMbsMutations() {
		return numNonPromoterMbsMutations;
	}

	public void setNumNonPromoterMbsMutations(int numNonPromoterMbsMutations) {
		this.numNonPromoterMbsMutations = numNonPromoterMbsMutations;
	}

	public int getNumPromoterMbsMutations() {
		return numPromoterMbsMutations;
	}

	public void setNumPromoterMbsMutations(int numPromoterMbsMutations) {
		this.numPromoterMbsMutations = numPromoterMbsMutations;
	}

	public int getNumPromoterMbsMutationsHighExp() {
		return numPromoterMbsMutationsHighExp;
	}

	public void setNumPromoterMbsMutationsHighExp(int numPromoterMbsMutationsHighExp) {
		this.numPromoterMbsMutationsHighExp = numPromoterMbsMutationsHighExp;
	}

	public int getNumPromoterMbsMutationsMidExp() {
		return numPromoterMbsMutationsMidExp;
	}

	public void setNumPromoterMbsMutationsMidExp(int numPromoterMbsMutationsMidExp) {
		this.numPromoterMbsMutationsMidExp = numPromoterMbsMutationsMidExp;
	}

	public int getNumPromoterMbsMutationsLowExp() {
		return numPromoterMbsMutationsLowExp;
	}

	public void setNumPromoterMbsMutationsLowExp(int numPromoterMbsMutationsLowExp) {
		this.numPromoterMbsMutationsLowExp = numPromoterMbsMutationsLowExp;
	}

	public int getNumPromoterMbsMutationsZeroExp() {
		return numPromoterMbsMutationsZeroExp;
	}

	public void setNumPromoterMbsMutationsZeroExp(int numPromoterMbsMutationsZeroExp) {
		this.numPromoterMbsMutationsZeroExp = numPromoterMbsMutationsZeroExp;
	}

	public int getNumNonCfsMutations() {
		return numNonCfsMutations;
	}

	public void setNumNonCfsMutations(int numNonCfsMutations) {
		this.numNonCfsMutations = numNonCfsMutations;
	}

	public int getNumCfsMutations() {
		return numCfsMutations;
	}

	public void setNumCfsMutations(int numCfsMutations) {
		this.numCfsMutations = numCfsMutations;
	}

	public int getNumNonCfsSbsMutations() {
		return numNonCfsSbsMutations;
	}

	public void setNumNonCfsSbsMutations(int numNonCfsSbsMutations) {
		this.numNonCfsSbsMutations = numNonCfsSbsMutations;
	}

	public int getNumCfsSbsMutations() {
		return numCfsSbsMutations;
	}

	public void setNumCfsSbsMutations(int numCfsSbsMutations) {
		this.numCfsSbsMutations = numCfsSbsMutations;
	}

	public int getNumNonCfsIndelMutations() {
		return numNonCfsIndelMutations;
	}

	public void setNumNonCfsIndelMutations(int numNonCfsIndelMutations) {
		this.numNonCfsIndelMutations = numNonCfsIndelMutations;
	}

	public int getNumCfsIndelMutations() {
		return numCfsIndelMutations;
	}

	public void setNumCfsIndelMutations(int numCfsIndelMutations) {
		this.numCfsIndelMutations = numCfsIndelMutations;
	}

	public int getNumNonCfsMbsMutations() {
		return numNonCfsMbsMutations;
	}

	public void setNumNonCfsMbsMutations(int numNonCfsMbsMutations) {
		this.numNonCfsMbsMutations = numNonCfsMbsMutations;
	}

	public int getNumCfsMbsMutations() {
		return numCfsMbsMutations;
	}

	public void setNumCfsMbsMutations(int numCfsMbsMutations) {
		this.numCfsMbsMutations = numCfsMbsMutations;
	}
	
	public String getLine() {
		List<Integer> lineElements = new ArrayList<>();
		
		lineElements.add(numSpecimens);
		lineElements.add(numSamples);
		lineElements.add(numSpecimensWithMutationData);
		lineElements.add(numSpecimensWithExpressionData);
		lineElements.add(numSpecimensWithBoth);
		
		lineElements.add(numMutations);
		
		lineElements.add(numNonPromoterMutations);
		lineElements.add(numPromoterMutations);
		
		lineElements.add(numNonPromoterSbsMutations);
		lineElements.add(numPromoterSbsMutations);
		lineElements.add(numPromoterSbsMutationsHighExp);
		lineElements.add(numPromoterSbsMutationsMidExp);
		lineElements.add(numPromoterSbsMutationsLowExp);
		lineElements.add(numPromoterSbsMutationsZeroExp);
		
		lineElements.add(numNonPromoterIndelMutations);
		lineElements.add(numPromoterIndelMutations);
		lineElements.add(numPromoterIndelMutationsHighExp);
		lineElements.add(numPromoterIndelMutationsMidExp);
		lineElements.add(numPromoterIndelMutationsLowExp);
		lineElements.add(numPromoterIndelMutationsZeroExp);
		
		lineElements.add(numNonPromoterMbsMutations);
		lineElements.add(numPromoterMbsMutations);
		lineElements.add(numPromoterMbsMutationsHighExp);
		lineElements.add(numPromoterMbsMutationsMidExp);
		lineElements.add(numPromoterMbsMutationsLowExp);
		lineElements.add(numPromoterMbsMutationsZeroExp);
		
		lineElements.add(numNonCfsMutations);
		lineElements.add(numCfsMutations);
		
		lineElements.add(numNonCfsSbsMutations);
		lineElements.add(numCfsSbsMutations);
		
		lineElements.add(numNonCfsIndelMutations);
		lineElements.add(numCfsIndelMutations);
		
		lineElements.add(numNonCfsMbsMutations);
		lineElements.add(numCfsMbsMutations);
		
		return StringUtils.join(lineElements, '\t');
	}
	
}
