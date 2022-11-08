package org.coh.mattiyeh.datacruncher.model;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import org.coh.mattiyeh.datacruncher.math.Operator;
import org.coh.mattiyeh.datacruncher.math.PercentileUtil;

public class Specimen {
	
	String specimenId;
	String donorId;
	String specimenType;
	String specimenSubtype;
	boolean tumorConfirmed;
	String donorTreatmentType;

	private Map<String, Sample> samples;

	public Specimen(String specimenId, String donorId, String specimenType, String specimenSubType, String donorTreatmentType) {
		super();
		this.specimenId = specimenId;
		this.donorId = donorId;
		this.specimenType = specimenType;
		this.specimenSubtype = specimenSubType;
		this.donorTreatmentType = donorTreatmentType;
		samples = new HashMap<>();
	}

	public String getSpecimenId() {
		return specimenId;
	}

	public String getDonorId() {
		return donorId;
	}

	public String getSpecimenType() {
		return specimenType;
	}

	public String getDonorTreatmentType() {
		return donorTreatmentType;
	}

	public boolean isTumorConfirmed() {
		return tumorConfirmed;
	}

	public Sample getSample(String sampleId) {
		return samples.get(sampleId);
	}

	public void addSample(Sample sample) {
		samples.put(sample.getSampleId(), sample);
	}

	private Sample getMutationSample() {
		for (Map.Entry<String, Sample> entry : samples.entrySet()) {
			Sample sample = entry.getValue();
			if (sample.hasMutationData()) {
				return sample;
			}
		}
		return null;
	}

	private Sample getExpressionSample() {
		for (Map.Entry<String, Sample> entry : samples.entrySet()) {
			Sample sample = entry.getValue();
			if (sample.hasExpressionData()) {
				return sample;
			}
		}
		return null;
	}

	/**
	 * It's okay if it's across two different SAMPLES
	 * 
	 * @return
	 */
	public boolean hasMutationAndExpressionData() {
		return hasMutationData() && hasExpressionData();
	}

	public boolean hasMutationData() {
		for (Map.Entry<String, Sample> entry : samples.entrySet()) {
			Sample sample = entry.getValue();
			if (sample.hasMutationData()) {
				return true;
			}
		}
		return false;
	}

	public boolean hasExpressionData() {
		for (Map.Entry<String, Sample> entry : samples.entrySet()) {
			Sample sample = entry.getValue();
			if (sample.hasExpressionData()) {
				return true;
			}
		}
		return false;
	}

	public int getNumSamples() {
		return samples.size();
	}

	public Map<String, Mutation> getMutations() {
		Map<String, Mutation> mutations = new HashMap<>();

		Sample mutationSample = getMutationSample();
		if (mutationSample != null) {
			mutations = mutationSample.getMutations();
		}

		return mutations;
	}

	public int getNumMutations() {
		return getMutations().size();
	}

	public Set<Mutation> getPromoterMutations() {
		Set<Mutation> promoterMutations = new TreeSet<>();

		Sample mutationSample = getMutationSample();
		if (mutationSample != null) {
			promoterMutations = mutationSample.getPromoterMutations();
		}

		return promoterMutations;
	}

	public int getNumPromoterMutations() {
		return getPromoterMutations().size();
	}
	
	public Set<Mutation> getNonPromoterMutations() {
		Set<Mutation> nonPromoterMutations = new TreeSet<>();

		Sample mutationSample = getMutationSample();
		if (mutationSample != null) {
			nonPromoterMutations = mutationSample.getNonPromoterMutations();
		}

		return nonPromoterMutations;
	}

	public int getNumNonPromoterMutations() {
		return getNonPromoterMutations().size();
	}

	public Set<Mutation> getCfsMutations() {
		Set<Mutation> cfsMutations = new TreeSet<>();

		Sample mutationSample = getMutationSample();
		if (mutationSample != null) {
			cfsMutations = mutationSample.getCfsMutations();
		}

		return cfsMutations;
	}

	public int getNumCfsMutations() {
		return getCfsMutations().size();
	}

	public Set<Mutation> getPromoterMutationsInExpressedGenes(int nthPercentile, Operator op) {
		Set<Mutation> promoterMutationsInExpressedGenes = new TreeSet<>();

		// Internal check to make sure at least one sample has mutation data and one
		// sample has expression data (it may be the same sample!)
		if (!hasMutationAndExpressionData()) {
			return promoterMutationsInExpressedGenes;
		}

		// First calculate the cutoff using the expression Sample
		double cutoff = 0;
		Sample expressionSample = getExpressionSample();
		if (expressionSample != null) {
			if (nthPercentile == 0) {
				cutoff = 0;
			} else {
				cutoff = PercentileUtil.calculateNthPercentile(expressionSample.getGeneNormExpressionLevelLogValues(),
						nthPercentile);
			}
		}

		// Now look for a sample with mutation data
		Sample mutationSample = getMutationSample();
		if (mutationSample != null) {
			promoterMutationsInExpressedGenes = mutationSample.getPromoterMutationsInExpressedGenes(cutoff, op,
					expressionSample);
		}

		return promoterMutationsInExpressedGenes;
	}

	public int getNumPromoterMutationsInExpressedGenes(int nthPercentile, Operator op) {
		return getPromoterMutationsInExpressedGenes(nthPercentile, op).size();
	}

	public boolean isPrimary() {
		return "Primary tumour".equals(specimenType);
	}

	public boolean isRecurrence() {
		return "Recurrent tumour".equals(specimenType);
	}

	public boolean isMetastasis() {
		return "Metastatic tumour".equals(specimenType);
	}

	public boolean isValidSpecimenSubtype() {
		return isSolidTissueSubtype() || "additional new primary".equals(specimenSubtype)
				|| "other".equals(specimenSubtype) || "lymph node".equals(specimenSubtype);
	}

	public boolean isSolidTissueSubtype() {
		return "solid tissue".equals(specimenSubtype);
	}

}
