package org.coh.mattiyeh.datacruncher.model;

import java.util.HashMap;
import java.util.List;
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
	boolean isTumorConfirmed;
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
		return isTumorConfirmed;
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
	
	public Set<Mutation> getMutations(MutationRange mutationRange, MutationType mutationType, Operator op, int nthPercentile, boolean useCutoff) {
		Set<Mutation> mutationsToReturn = new TreeSet<>();
		Sample expressionSample = null;
		double cutoff = 0;
		
		// Internal check to make sure at least one sample has mutation data
		if (!hasMutationData()) {
			return mutationsToReturn;
		}
		
		if (useCutoff) {
			
			// Internal check to make sure at least one sample has expression data (it
			// doesn't have to be the same sample as mutation sample above!)
			if (!hasExpressionData()) {
				return mutationsToReturn;
			}

			// First calculate the cutoff using the expression Sample
			expressionSample = getExpressionSample();
			
			if (expressionSample != null) {
				if (nthPercentile == 0) {
					cutoff = 0;
				} else {
					List<Double> geneNormExpressionLevelLogValues = expressionSample
							.getGeneNormExpressionLevelLogValues();
					cutoff = PercentileUtil.calculateNthPercentile(geneNormExpressionLevelLogValues, nthPercentile);
				}
			} else {
				// TODO: throw exception if null?
			}
			
		}
		
		// Now look for a sample with mutation data
		Sample mutationSample = getMutationSample();
		
		if (mutationSample != null) {
			mutationsToReturn = mutationSample.getMutations(mutationRange, mutationType, op, cutoff, expressionSample);
		}
		
		return mutationsToReturn;
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
