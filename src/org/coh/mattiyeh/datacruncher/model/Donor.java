package org.coh.mattiyeh.datacruncher.model;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

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

	private Specimen pickBestSpecimen(Specimen oldSpecimen, Specimen newSpecimen) {
		Specimen specimenToReturn = oldSpecimen;

		// Reasons to get a new specimen:
		// 1. New one is primary when old one is not
		// 2. Both are primary but new one is a BETTER primary (eg. solid tissue)
		// 3. New one is recurrence when old one is met

		if (newSpecimen.isPrimary()) {
			if (!oldSpecimen.isPrimary()) {
				// Scenario #1
				specimenToReturn = newSpecimen;
			} else {
				if (newSpecimen.isSolidTissueSubtype() && !oldSpecimen.isSolidTissueSubtype()) {
					// Scenario #2
					specimenToReturn = newSpecimen;
				}
			}
		} else if (newSpecimen.isRecurrence() && oldSpecimen.isMetastasis()) {
			// Scenario #3
			specimenToReturn = newSpecimen;
		}
		return specimenToReturn;
	}
	
	/**
	 * This makes sure that at least one specimen has mutation data.
	 * 
	 * @return
	 */
	public boolean hasMutationData() {
		for (Map.Entry<String, Specimen> entry : specimens.entrySet()) {
			Specimen specimen = entry.getValue();
			if (specimen.hasMutationData()) {
				return true;
			}
		}
		// If we reach here, it means we never returned true in the loop meaning NONE of
		// the specimens have mutation data.
		return false;
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

	public Set<Mutation> getIndelMutations(MutationRange mutationRange) {
		Set<Mutation> mutations = getMutations(mutationRange, MutationType.INSERTION);
		mutations.addAll(getMutations(mutationRange, MutationType.DELETION));
		return mutations;
	}
	
	public Set<Mutation> getIndelMutationsInExpressedGenes(MutationRange mutationRange, Operator op, int nthPercentile) {
		Set<Mutation> mutations = getMutations(mutationRange, MutationType.INSERTION, op, nthPercentile, true);
		mutations.addAll(getMutations(mutationRange, MutationType.DELETION, op, nthPercentile, true));
		return mutations;
	}
	
	public Set<Mutation> getIndelMutationsInExpressedGenes(MutationRange mutationRange, Operator lowOp, int lowNthPercentile, Operator highOp, int highNthPercentile) {
		Set<Mutation> mutations = getMutationsInExpressedGenes(mutationRange, MutationType.INSERTION, lowOp, lowNthPercentile, highOp, highNthPercentile);
		mutations.addAll(getMutationsInExpressedGenes(mutationRange, MutationType.DELETION, lowOp, lowNthPercentile, highOp, highNthPercentile));
		return mutations;
	}

	public Set<Mutation> getMutations(MutationRange mutationRange, MutationType mutationType) {
		return getMutations(mutationRange, mutationType, Operator.GREATERTHANOREQUAL, 0, false);
	}

	public Set<Mutation> getMutationsInExpressedGenes(MutationRange mutationRange, MutationType mutationType, Operator op, int nthPercentile) {
		return getMutations(mutationRange, mutationType, op, nthPercentile, true);
	}

	/**
	 * @param mutationRange
	 * @param mutationType
	 * @param lowOp
	 * @param lowNthPercentile
	 * @param highOp
	 * @param highNthPercentile
	 * @return
	 */
	public Set<Mutation> getMutationsInExpressedGenes(MutationRange mutationRange, MutationType mutationType, Operator lowOp, int lowNthPercentile, Operator highOp, int highNthPercentile) {
		// Gets all mutations greater than 25%
		Set<Mutation> lowMutations = getMutations(mutationRange, mutationType, lowOp, lowNthPercentile, true);
		
		// Gets all mutations less than 75%
		Set<Mutation> highMutations = getMutations(mutationRange, mutationType, highOp, highNthPercentile, true);
		
		// Intersect both sets to keep only those that are between 25% and 75%
		lowMutations.retainAll(highMutations);
		
		return lowMutations;
	}

	/**
	 * @param mutationRange
	 * @param mutationType
	 * @param op
	 * @param nthPercentile
	 * @param useCutoff
	 * @return
	 */
	public Set<Mutation> getMutations(MutationRange mutationRange, MutationType mutationType, Operator op,
			int nthPercentile, boolean useCutoff) {
		Set<Mutation> mutations = new TreeSet<>();

		Specimen specimen = findSpecimenWithMutationAndExpressionData();
		if (specimen != null) {
			mutations = specimen.getMutations(mutationRange, mutationType, op, nthPercentile, useCutoff);
		}

		return mutations;
	}

	public boolean containsSpecimen(String specimenId) {
		return specimens.containsKey(specimenId);
	}

}
