package org.coh.mattiyeh.datacruncher.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

import org.apache.commons.lang3.StringUtils;
import org.coh.mattiyeh.datacruncher.Constants;
import org.coh.mattiyeh.datacruncher.genome.GenomeUtils;
import org.coh.mattiyeh.datacruncher.genome.TriSeq;

public class Mutation implements Serializable, Comparable<Mutation> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 277334304446372552L;

	private String mutationId;
	private String donorId;
	private String specimenId;
	private String sampleId;
	private String matchedSampleId;
	private MutationType type;
	private String chr;
	private int start;
	private int end;
	private String refBase;
	private String mutBase;
	private int totalReadCount;
	private int mutantAlleleReadCount;
	private String sequencingStrategy;
	private List<MutationEffect> mutationEffects = new ArrayList<>();

	private TriSeq triSeq;

	private boolean isInPromoterRegion;
	private boolean isInCfsRegion;

	public Mutation(String mutationId, String donorId, String specimenId, String sampleId, String matchedSampleId,
			String type, String chr, int start, int end, String refBase, String mutBase, int totalReadCount,
			int mutantAlleleReadCount, String sequencingStrategy) {

		this.mutationId = mutationId;
		this.donorId = donorId;
		this.specimenId = specimenId;
		this.sampleId = sampleId;
		this.matchedSampleId = matchedSampleId;

		this.type = switch (type) {
		case Constants.SBS			-> MutationType.SBS;
		case Constants.INSERTION	-> MutationType.INSERTION;
		case Constants.DELETION		-> MutationType.DELETION;
		case Constants.MBS			-> MutationType.MBS;
		default -> throw new IllegalStateException("Invalid type: " + type);
		};

		this.chr = chr;
		this.start = start;
		this.end = end;
		this.refBase = refBase;
		this.mutBase = mutBase;
		this.totalReadCount = totalReadCount;
		this.mutantAlleleReadCount = mutantAlleleReadCount;
		this.sequencingStrategy = sequencingStrategy;

		isInPromoterRegion = false;
		isInCfsRegion = false;

	}

	public String getMutationId() {
		return mutationId;
	}

	public String getDonorId() {
		return donorId;
	}

	public String getSpecimenId() {
		return specimenId;
	}

	public String getSampleId() {
		return sampleId;
	}

	public String getMatchedSampleId() {
		return matchedSampleId;
	}

	public MutationType getType() {
		return type;
	}

	public TriSeq getTriSeq() {
		return triSeq;
	}

	public void setTriSeq(TriSeq triSeq) {
		this.triSeq = triSeq;
	}

	public String getChr() {
		return chr;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public String getRefBase() {
		return refBase;
	}

	public String getMutBase() {
		return mutBase;
	}

	public int getTotalReadCount() {
		return totalReadCount;
	}

	public int getMutantAlleleReadCount() {
		return mutantAlleleReadCount;
	}

	public String getSequencingStrategy() {
		return sequencingStrategy;
	}

	public List<MutationEffect> getMutationEffects() {
		return mutationEffects;
	}

	public boolean isInPromoterRegion() {
		return isInPromoterRegion;
	}

	public void setIsInPromoterRegion(boolean isInPromoterRegion) {
		this.isInPromoterRegion = isInPromoterRegion;
	}

	public boolean isInCfsRegion() {
		return isInCfsRegion;
	}

	public void setIsInCfsRegion(boolean isInCfsRegion) {
		this.isInCfsRegion = isInCfsRegion;
	}

	public void addMutationEffect(MutationEffect mutationEffect) {
		mutationEffects.add(mutationEffect);
	}

	public List<String> getGeneIdsAffected() {
		List<String> genesAffected = new ArrayList<>();

		getProteinCodingMutationEffects()
				.forEach(mutationEffect -> genesAffected.add(mutationEffect.getGeneAffected()));

		return genesAffected;
	}
	
	private List<MutationEffect> getProteinCodingMutationEffects() {
		List<MutationEffect> proteinCodingMutationEffects = new ArrayList<>();
		for (MutationEffect mutationEffect : mutationEffects) {
			if (mutationEffect.isProteinCoding()) {
				proteinCodingMutationEffects.add(mutationEffect);
			}
		}
		return proteinCodingMutationEffects;
	}

	public boolean affectsProteinCodingGene() {
		for (MutationEffect mutationEffect : mutationEffects) {
			if (mutationEffect.isProteinCoding()) {
				return true;
			}
		}
		return false;
	}

	public boolean affectsGene(String geneId) {
		if (StringUtils.isBlank(geneId)) {
			return false;
		}
		
		for (MutationEffect mutationEffect : mutationEffects) {
			if (mutationEffect.getGeneAffected().equals(geneId)) {
				return true;
			}
		}
		return false;
	}

	public int getNumGeneIdsAffected() {
		return getGeneIdsAffected().size();
	}
	
	public String getFirstRawLine() {
		List<MutationEffect> proteinCodingMutationEffects = getProteinCodingMutationEffects();
		if (proteinCodingMutationEffects.isEmpty()) {
			return getMutationEffects().get(0).getRawLine(); 
		}
		return proteinCodingMutationEffects.get(0).getRawLine();
	}

	public String getTriSeqWithMut() {
		if (triSeq == null) {
			return StringUtils.EMPTY;
		}

		String change = "T:" + triSeq.getPreBase() + "[" + triSeq.getRefBase() + ">" + mutBase + "]"
				+ triSeq.getPostBase();

		return change.toUpperCase();
	}

	public String getTriSeqWithMutForSigs() {
		if (triSeq == null) {
			return StringUtils.EMPTY;
		}

		// If it's already in good form, then return the original function
		if ("C".equals(refBase) || "T".equals(refBase)) {
			return getTriSeqWithMut();
		}

		// Prebase and postbase are SWITCHED since it's the complement strand and it's
		// read backwards
		// E.g. C[G>A]T becomes A[C>T]G
		String preBaseToWrite = triSeq.getPreBase();
		String refBaseToWrite = refBase;
		String mutBaseToWrite = mutBase;
		String postBaseToWrite = triSeq.getPostBase();

		// Need to get reverse complement in order to turn A> and G> into C> and T>

		preBaseToWrite = GenomeUtils.getReverseComplement(preBaseToWrite);
		refBaseToWrite = GenomeUtils.getReverseComplement(refBaseToWrite);
		mutBaseToWrite = GenomeUtils.getReverseComplement(mutBaseToWrite);
		postBaseToWrite = GenomeUtils.getReverseComplement(postBaseToWrite);

		// NOT A TYPO. Read above.
		String change = "T:" + postBaseToWrite + "[" + refBaseToWrite + ">" + mutBaseToWrite + "]" + preBaseToWrite;

		return change.toUpperCase();
	}

	@Override
	public int hashCode() {
		return Objects.hash(mutationId);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Mutation other = (Mutation) obj;
		return Objects.equals(mutationId, other.mutationId);
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		builder.append("Mutation [mutationId=");
		builder.append(mutationId);
		builder.append(", chr=");
		builder.append(chr);
		builder.append(", start=");
		builder.append(start);
		builder.append(", end=");
		builder.append(end);
		builder.append(", refBase=");
		builder.append(refBase);
		builder.append(", mutBase=");
		builder.append(mutBase);
		builder.append("]");
		return builder.toString();
	}

	@Override
	public int compareTo(Mutation o) {
		return this.mutationId.compareTo(o.getMutationId());
	}

}
