package org.coh.mattiyeh.datacruncher.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

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

	public boolean hasPromoterMutation() {
		return getNumPromoterMutations() > 0;
	}

	public boolean hasPromoterMutationAndExpressionData() {
		return hasMutationAndExpressionData() && hasPromoterMutation();
	}

	public Map<String, Mutation> getMutations() {
		return mutations;
	}

	public Mutation getMutation(String mutationId) {
		return mutations.get(mutationId);
	}

	public int getNumMutations() {
		return mutations.size();
	}

	public List<Mutation> getPromoterMutations() {
		List<Mutation> promoterMutations = new ArrayList<>();

		mutations.forEach((mutationId, mutation) -> {
			if (mutation.isInPromoterRegion()) {
				promoterMutations.add(mutation);
			}
		});

		return promoterMutations;
	}
	
	public int getNumPromoterMutations() {
		return getPromoterMutations().size();
	}

	public List<Mutation> getNonPromoterMutations() {
		List<Mutation> nonPromoterMutations = new ArrayList<>();

		mutations.forEach((mutationId, mutation) -> {
			if (!mutation.isInPromoterRegion()) {
				nonPromoterMutations.add(mutation);
			}
		});

		return nonPromoterMutations;
	}

	public int getNumNonPromoterMutations() {
		return getNonPromoterMutations().size();
	}

	public List<Mutation> getCfsMutations() {
		List<Mutation> cfsMutations = new ArrayList<>();

		mutations.forEach((mutationId, mutation) -> {
			if (mutation.isInCfsRegion()) {
				cfsMutations.add(mutation);
			}
		});

		return cfsMutations;
	}

	public int getNumCfsMutations() {
		return getCfsMutations().size();
	}

	public Set<Mutation> getPromoterMutationsInExpressedGenes(double cutoff, Operator op, Sample expressionSample) {
		Set<Mutation> promoterMutationsInExpressedGenes = new HashSet<>();

		// All gene IDs in this sample we have expression data for
		Set<String> geneIds = expressionSample.getGeneNormExpressionLevels().keySet();
		for (String geneId : geneIds) {
			// Is gene highly expressed? ... or whatever the operator calls for
			if (op.apply(expressionSample.getGeneNormExpressionLevelLog(geneId), cutoff)) {

				// Is there a promoter mutation for this gene?
				for (Mutation mutation : getMutationsForGene(geneId)) {
					if (mutation.isInPromoterRegion()) {
						promoterMutationsInExpressedGenes.add(mutation);
					}
				}
			}
		}

		return promoterMutationsInExpressedGenes;
	}

	public int getNumPromoterMutationsInExpressedGenes(double cutoff, Operator op, Sample expressionSample) {
		return getPromoterMutationsInExpressedGenes(cutoff, op, expressionSample).size();
	}

	public void addGeneNormExpressionLevel(String geneId, float normalizedReadCount) {
		geneNormExpressionLevels.put(geneId, normalizedReadCount);
	}

	public Map<String, Float> getGeneNormExpressionLevels() {
		return geneNormExpressionLevels;
	}

	public Double getGeneNormExpressionLevelLog(String geneId) {
		return Functions.logTransform(MapUtils.getFloat(geneNormExpressionLevels, geneId, 0f));
	}

	public List<Double> getGeneNormExpressionLevelLogValues() {
		List<Double> geneNormExpressionLevelLogValues = new ArrayList<>();
		geneNormExpressionLevels.forEach(
				(geneId, value) -> geneNormExpressionLevelLogValues.add(getGeneNormExpressionLevelLog(geneId)));
		return geneNormExpressionLevelLogValues;
	}

	private Set<Mutation> getMutationsForGene(String geneId) {
		Set<Mutation> mutationsForGene = new HashSet<>();

		mutations.forEach((mutationId, mutation) -> {
			if (mutation.affectsGene(geneId)) {
				mutationsForGene.add(mutation);
			}
		});

		return mutationsForGene;
	}

}
