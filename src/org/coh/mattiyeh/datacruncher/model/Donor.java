package org.coh.mattiyeh.datacruncher.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.coh.mattiyeh.datacruncher.math.Operator;

public class Donor {

	String donorId;
	int survivalTime;

	private Map<String, Specimen> specimens;

	public Donor(String donorId, int survivalTime) {
		super();
		this.donorId = donorId;
		this.survivalTime = survivalTime;

		specimens = new HashMap<>();
	}

	public String getDonorId() {
		return donorId;
	}
	
	public Specimen getSpecimen(String specimenId) {
		return specimens.get(specimenId);
	}

	public Map<String, Specimen> getSpecimens() {
		return specimens;
	}

	public void addSpecimen(Specimen specimen) {
		specimens.put(specimen.getSpecimenId(), specimen);
	}

	/**
	 * @return
	 */
	private Specimen findSpecimenWithMutationData() {
		Specimen specimenToReturn = null;

		// Note: all normal, cell line, and xenografts were already filtered out during
		// read-in (See DataCruncher#validSpecimenType)

		// See https://docs.icgc.org/submission/guide/sample-type-specification/ for
		// specimen types

		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen newSpecimen = entry.getValue();

			if (!newSpecimen.isValidSpecimenSubtype() || !newSpecimen.hasMutationData()) {
				continue;
			}

			// At this point, all specimens are good subtype (ie. solid tissue, other, addl
			// new primary)

			if (specimenToReturn == null) {
				// First candidate
				specimenToReturn = newSpecimen;
			} else {

				// Reasons to get a new specimen:
				// 1. New one is primary when old one is not
				// 2. Both are primary but new one is a BETTER primary (eg. solid tissue)
				// 3. New one is recurrence when old one is met

				specimenToReturn = pickBestSpecimen(specimenToReturn, newSpecimen);
			}
		}

		return specimenToReturn;
	}

	private Specimen findSpecimenWithMutationAndExpressionData() {
		Specimen specimenToReturn = null;

		// Note: all normal, cell line, and xenografts were already filtered out during
		// read-in (See DataCruncher#validSpecimenType)

		// See https://docs.icgc.org/submission/guide/sample-type-specification/ for
		// specimen types

		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen newSpecimen = entry.getValue();

			if (!newSpecimen.isValidSpecimenSubtype() || !newSpecimen.hasMutationAndExpressionData()) {
				continue;
			}

			// At this point, all specimens are good subtype (ie. solid tissue, other, addl
			// new primary)

			if (specimenToReturn == null) {
				// First candidate
				specimenToReturn = newSpecimen;
			} else {

				specimenToReturn = pickBestSpecimen(specimenToReturn, newSpecimen);
			}
		}

		// If it can't find the perfect specimen, then return a specimen with just
		// mutation data
		// This fixes a bug where two specimens have mutation data and different ones
		// get picked depending on context
		if (specimenToReturn == null) {
			return findSpecimenWithMutationData();
		}

		return specimenToReturn;
	}

	private Specimen pickBestSpecimen(Specimen specimenToReturn, Specimen newSpecimen) {
		// Reasons to get a new specimen:
		// 1. New one is primary when old one is not
		// 2. Both are primary but new one is a BETTER primary (eg. solid tissue)
		// 3. New one is recurrence when old one is met

		if (newSpecimen.isPrimary()) {
			if (!specimenToReturn.isPrimary()) {
				// Scenario #1
				specimenToReturn = newSpecimen;
			} else {
				if (newSpecimen.isSolidTissueSubtype() && !specimenToReturn.isSolidTissueSubtype()) {
					// Scenario #2
					specimenToReturn = newSpecimen;
				}
			}
		} else if (newSpecimen.isRecurrence() && specimenToReturn.isMetastasis()) {
			// Scenario #3
			specimenToReturn = newSpecimen;
		}
		return specimenToReturn;
	}

	/**
	 * This makes sure that at least one specimen has BOTH mutation and expression
	 * data. Note: the mutation and expression data don't have to come from the same
	 * SAMPLE, but they should be from same SPECIMEN
	 * 
	 * @return
	 */
	public boolean hasMutationAndExpressionData() {
		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen specimen = entry.getValue();
			if (specimen.hasMutationAndExpressionData()) {
				return true;
			}
		}
		// If we reach here, it means we never returned true in the loop meaning NONE of
		// the specimens had both mutation and expression data.
		return false;
	}

	public int getNumSpecimens() {
		return specimens.size();
	}

	public int getNumSpecimensWithMutationData() {
		int total = 0;
		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen specimen = entry.getValue();
			if (specimen.hasMutationData()) {
				total++;
			}
		}
		return total;
	}

	public int getNumSpecimensWithExpressionData() {
		int total = 0;
		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen specimen = entry.getValue();
			if (specimen.hasExpressionData()) {
				total++;
			}
		}
		return total;
	}

	public int getNumSpecimensWithBoth() {
		int total = 0;
		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen specimen = entry.getValue();
			if (specimen.hasMutationAndExpressionData()) {
				total++;
			}
		}
		return total;
	}

	public int getNumSamples() {
		int total = 0;
		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen specimen = entry.getValue();
			total += specimen.getNumSamples();
		}
		return total;
	}

	public Map<String, Mutation> getMutations() {
		Map<String, Mutation> mutations = new HashMap<>();

		Specimen specimen = findSpecimenWithMutationAndExpressionData();
		if (specimen != null) {
			mutations = specimen.getMutations();
		}

		return mutations;
	}

	public int getNumMutations() {
		return getMutations().size();
	}

	public List<Mutation> getPromoterMutations() {
		List<Mutation> promoterMutations = new ArrayList<>();

		Specimen specimen = findSpecimenWithMutationAndExpressionData();
		if (specimen != null) {
			promoterMutations = specimen.getPromoterMutations();
		}

		return promoterMutations;
	}

	public int getNumPromoterMutations() {
		return getPromoterMutations().size();
	}
	
	public List<Mutation> getNonPromoterMutations() {
		List<Mutation> nonPromoterMutations = new ArrayList<>();

		Specimen specimen = findSpecimenWithMutationAndExpressionData();
		if (specimen != null) {
			nonPromoterMutations = specimen.getNonPromoterMutations();
		}

		return nonPromoterMutations;
	}

	public int getNumNonPromoterMutations() {
		return getNonPromoterMutations().size();
	}
	
	public List<Mutation> getCfsMutations() {
		List<Mutation> cfsMutations = new ArrayList<>();

		Specimen specimen = findSpecimenWithMutationAndExpressionData();
		if (specimen != null) {
			cfsMutations = specimen.getCfsMutations();
		}

		return cfsMutations;
	}

	public int getNumCfsMutations() {
		return getCfsMutations().size();
	}

	public Set<Mutation> getPromoterMutationsInExpressedGenes(int nthPercentile, Operator op) {
		Set<Mutation> promoterMutationsInExpressedGenes = new HashSet<>();

		Specimen specimen = findSpecimenWithMutationAndExpressionData();
		if (specimen != null) {
			promoterMutationsInExpressedGenes = specimen.getPromoterMutationsInExpressedGenes(nthPercentile, op);
		}

		return promoterMutationsInExpressedGenes;

	}

	public int getNumPromoterMutationsInHighlyExpressedGenes(int nthPercentile, Operator op) {
		return getPromoterMutationsInExpressedGenes(nthPercentile, op).size();
	}

	public boolean containsSpecimen(String specimenId) {
		return specimens.containsKey(specimenId);
	}

}
