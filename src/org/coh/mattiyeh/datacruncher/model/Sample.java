package org.coh.mattiyeh.datacruncher.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.apache.commons.collections4.MapUtils;
import org.coh.mattiyeh.datacruncher.math.Functions;
import org.coh.mattiyeh.datacruncher.math.Operator;

public class Sample {

	String sampleId;
	String specimenId;
	String donorId;

	private Map<String, Mutation> mutations;
	private Map<String, Float> geneNormExpressionLevels;

	public Sample(String sampleId, String specimenId, String donorId) {
		this(sampleId, specimenId, donorId, new HashMap<>(), new HashMap<>());
	}

	public Sample(String sampleId, String specimenId, String donorId, Map<String, Mutation> mutations, Map<String, Float> geneNormExpressionLevels) {
		super();
		this.sampleId = sampleId;
		this.specimenId = specimenId;
		this.donorId = donorId;

		this.mutations = mutations;
		this.geneNormExpressionLevels = geneNormExpressionLevels;
	}

	public String getSampleId() {
		return sampleId;
	}

	public String getSpecimenId() {
		return specimenId;
	}

	public void addMutation(Mutation mutation) {
		mutations.put(mutation.getMutationId(), mutation);
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

	public Map<String, Mutation> getMutations() {
		return mutations;
	}

	public Mutation getMutation(String mutationId) {
		return mutations.get(mutationId);
	}

	/**
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

			// Check to see if mutation is located in promoter/cfs range, outside
			// promoter/cfs, or we don't care (ie. NONE)
			// AND
			// If this mutation matches the type we're looking for (eg. SBS) or we're
			// looking for all types

			if (checkRange(mutationRange, mutation) && checkType(mutationType, mutation)) {

				// If expressionSample isn't null, then we have a request for mutations that
				// meet an expression criteria
				if (expressionSample == null) {
					mutationsToReturn.add(mutation);
				} else {
					mutation.getGeneIdsAffected().forEach(geneId -> {
						// Is gene highly expressed? ... or whatever the operator calls for
						if (op.apply(expressionSample.getGeneNormExpressionLevelLog(geneId), expressionCutoff)) {
							mutationsToReturn.add(mutation);
						}
					});
				}

			}

		});

		return mutationsToReturn;
	}

	private boolean checkRange(MutationRange mutationRange, Mutation mutation) {
		return (MutationRange.PROMOTER.equals(mutationRange) && mutation.isInPromoterRegion())
				|| (MutationRange.NONPROMOTER.equals(mutationRange) && !mutation.isInPromoterRegion())
				|| (MutationRange.CFS.equals(mutationRange) && mutation.isInCfsRegion())
				|| (MutationRange.NONCFS.equals(mutationRange) && !mutation.isInCfsRegion())
				|| (MutationRange.NONE.equals(mutationRange));
	}

	private boolean checkType(MutationType mutationType, Mutation mutation) {
		return mutationType.equals(mutation.getType()) || MutationType.ALL.equals(mutationType);
	}

	public void addGeneNormExpressionLevel(String geneId, float normalizedReadCount) {
		geneNormExpressionLevels.put(geneId, normalizedReadCount);
	}

	public Map<String, Float> getGeneNormExpressionLevels() {
		return geneNormExpressionLevels;
	}

	public Double getGeneNormExpressionLevelLog(String geneId) {
		return Functions.safeLogTransform(MapUtils.getFloat(geneNormExpressionLevels, geneId, 0f));
	}

	public List<Double> getGeneNormExpressionLevelLogValues() {
		List<Double> geneNormExpressionLevelLogValues = new ArrayList<>();
		geneNormExpressionLevels.forEach(
				(geneId, value) -> geneNormExpressionLevelLogValues.add(getGeneNormExpressionLevelLog(geneId)));
		return geneNormExpressionLevelLogValues;
	}

}
