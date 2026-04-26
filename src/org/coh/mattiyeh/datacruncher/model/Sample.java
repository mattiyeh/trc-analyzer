package org.coh.mattiyeh.datacruncher.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.coh.mattiyeh.datacruncher.math.Functions;
import org.coh.mattiyeh.datacruncher.math.Operator;

public class Sample {

	String sampleId;
	String specimenId;
	String donorId;

	private Map<String, Mutation> mutations = new HashMap<>();
	private Map<String, Float> geneNormExpressionLevels = new HashMap<>();

	public Sample(String sampleId, String specimenId, String donorId) {
		super();
		this.sampleId = sampleId;
		this.specimenId = specimenId;
		this.donorId = donorId;
	}

	public String getSampleId() {
		return sampleId;
	}

	public String getSpecimenId() {
		return specimenId;
	}

	public Map<String, Mutation> getMutations() {
		return mutations;
	}

	/**
	 * Returns the subset of this sample's mutations matching the requested
	 * range / type / expression criteria.
	 *
	 * <h3>Multi-gene attribution semantics</h3>
	 * If {@code expressionSample} is non-null, a mutation that overlaps the
	 * promoters of MULTIPLE protein-coding genes is included if <em>any</em>
	 * of those genes satisfies the expression cutoff under the given operator.
	 * The mutation is only added once (TreeSet deduplicates), but the same
	 * mutation can be returned by multiple separate calls with different
	 * operators (e.g. {@code GREATERTHANOREQUAL p75} and
	 * {@code LESSTHANOREQUAL p25}) when its overlapping genes fall on opposite
	 * sides of the cutoff.
	 *
	 * <p>Practical consequence: per-bin output files
	 * (*_high.tsv / *_mid.tsv / *_low.tsv / *_zero.tsv) DO NOT have row counts
	 * that sum to the unfiltered *_all.tsv count. A single mutation can appear
	 * in two bin files if it overlaps two genes whose expressions fall in
	 * different bins. Manuscript Methods should describe this explicitly when
	 * citing per-bin mutation counts.
	 *
	 * @param mutationRange
	 * @param mutationType
	 * @param op
	 * @param expressionCutoff
	 * @param expressionSample
	 * @return
	 */
	public Set<Mutation> getMutations(MutationRange mutationRange, MutationType mutationType, Operator op,
			double expressionCutoff, Sample expressionSample) {
		Set<Mutation> mutationsToReturn = new TreeSet<>();
	
		mutations.forEach((mutationId, mutation) -> {
			
			// Check to see if mutation has at least one mutationEffect that is protein
			// coding. If not, we don't count it and move on to next mutation.
			// The only time we want to count these mutations is if we're looking for
			// unfiltered mutations (hence the second part of the AND statement)
			if (!mutation.affectsProteinCodingGene() && !MutationRange.UNFILTERED.equals(mutationRange)) {
				return;
			}
			
			// Check to see if mutation is located in promoter/cfs range, outside
			// promoter/cfs, or we don't care (ie. NONE)
			// AND
			// If this mutation matches the type we're looking for (eg. SBS) or we're
			// looking for all types
	
			if (isValidRange(mutationRange, mutation) && isValidType(mutationType, mutation)) {
	
				// If expressionSample is null, then that means we have a request for mutations
				// that doesn't need expression data
				
				// If expressionSample isn't null, then we have a request for mutations that
				// meet an expression criteria
				
				if (expressionSample == null) {
					mutationsToReturn.add(mutation);
				} else {
					mutation.getGeneIdsAffected().forEach(geneId -> {
						
						// Make sure we have expression data for this gene:
						if (expressionSample.hasExpressionForGene(geneId)) {
							
							// Is gene highly expressed? ... or whatever the operator calls for
							if (op.apply(expressionSample.getGeneNormExpressionLevelLog(geneId), expressionCutoff)) {
								mutationsToReturn.add(mutation);
							}
						}
					});
				}
			}
		});
	
		return mutationsToReturn;
	}

	public void addMutation(Mutation mutation) {
		mutations.put(mutation.getMutationId(), mutation);
	}

	public Mutation getMutation(String mutationId) {
		return mutations.get(mutationId);
	}

	public boolean hasMutationData() {
		return !mutations.isEmpty();
	}

	public boolean hasExpressionData() {
		return !geneNormExpressionLevels.isEmpty();
	}

	public boolean hasMutationAndExpressionData() {
		return hasMutationData() && hasExpressionData();
	}

	private boolean hasExpressionForGene(String geneId) {
		return geneNormExpressionLevels.containsKey(geneId);
	}

	/**
	 * Check to see if mutation is in the range we care about (eg. promoter, cfs, either, none)
	 * 
	 * @param mutationRange
	 * @param mutation
	 * @return
	 */
	private boolean isValidRange(MutationRange mutationRange, Mutation mutation) {
		return (MutationRange.PROMOTER.equals(mutationRange) && mutation.isInPromoterRegion())
				|| (MutationRange.NONPROMOTER.equals(mutationRange) && !mutation.isInPromoterRegion())
				|| (MutationRange.CFS.equals(mutationRange) && mutation.isInCfsRegion())
				|| (MutationRange.NONCFS.equals(mutationRange) && !mutation.isInCfsRegion())
				|| (MutationRange.UNFILTERED.equals(mutationRange))
				|| (MutationRange.NONE.equals(mutationRange));
	}

	/**
	 * Check to see if mutation is the type we care about (ie. SBS, INSERTION, DELETION, MBS)
	 * 
	 * @param mutationType
	 * @param mutation
	 * @return
	 */
	private boolean isValidType(MutationType mutationType, Mutation mutation) {
		return mutationType.equals(mutation.getType()) || MutationType.ALL.equals(mutationType);
	}

	public void addGeneNormExpressionLevel(String geneId, float normalizedReadCount) {
		geneNormExpressionLevels.put(geneId, normalizedReadCount);
	}

	public Map<String, Float> getGeneNormExpressionLevels() {
		return geneNormExpressionLevels;
	}

	public Double getGeneNormExpressionLevelLog(String geneId) {
		return Functions.safeLogTransform(geneNormExpressionLevels.get(geneId));
	}

	public List<Double> getGeneNormExpressionLevelLogValues() {
		List<Double> geneNormExpressionLevelLogValues = new ArrayList<>();
		geneNormExpressionLevels.forEach(
				(geneId, value) -> geneNormExpressionLevelLogValues.add(getGeneNormExpressionLevelLog(geneId)));
		return geneNormExpressionLevelLogValues;
	}

}
