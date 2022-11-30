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
	private List<MutationEffect> mutationEffects;

	private TriSeq triSeq;
	
	private boolean inPromoterRegion;
	private boolean inCfsRegion;

	private List<String> rawLines;

	public Mutation(String mutationId, String donorId, String specimenId, String sampleId, String matchedSampleId,
			String type, String chr, int start, int end, String refBase, String mutBase, int totalReadCount,
			int mutantAlleleReadCount, String sequencingStrategy, String consequenceType, String geneAffected) {

		this.mutationId = mutationId;
		this.donorId = donorId;
		this.specimenId = specimenId;
		this.sampleId = sampleId;
		this.matchedSampleId = matchedSampleId;

		if (Constants.SBS.equals(type)) {
			this.type = MutationType.SBS;
		} else if (Constants.INSERTION.equals(type)) {
			this.type = MutationType.INSERTION;
		} else if (Constants.DELETION.equals(type)) {
			this.type = MutationType.DELETION;
		} else if (Constants.MBS.equals(type)) {
			this.type = MutationType.MBS;
		}

		this.chr = chr;
		this.start = start;
		this.end = end;
		this.refBase = refBase;
		this.mutBase = mutBase;
		this.totalReadCount = totalReadCount;
		this.mutantAlleleReadCount = mutantAlleleReadCount;
		this.sequencingStrategy = sequencingStrategy;

		mutationEffects = new ArrayList<>();
		mutationEffects.add(new MutationEffect(consequenceType, geneAffected));

		inPromoterRegion = false;
		inCfsRegion = false;

		this.rawLines = new ArrayList<>();
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
		return inPromoterRegion;
	}

	public void setInPromoterRegion(boolean inPromoterRegion) {
		this.inPromoterRegion = inPromoterRegion;
	}

	public boolean isInCfsRegion() {
		return inCfsRegion;
	}

	public void setInCfsRegion(boolean inCfsRegion) {
		this.inCfsRegion = inCfsRegion;
	}

	public void addMutationEffect(String consequenceType, String geneAffected) {
		mutationEffects.add(new MutationEffect(consequenceType, geneAffected));
	}

	public List<String> getGeneIdsAffected() {
		List<String> genesAffected = new ArrayList<>();
		mutationEffects.forEach(mutationEffect -> genesAffected.add(mutationEffect.getGeneAffected()));
		return genesAffected;
	}
	
	public int getNumGeneIdsAffected() {
		return getGeneIdsAffected().size();
	}

	public boolean affectsGene(String geneId) {
		for (MutationEffect mutationEffect : mutationEffects) {
			if (mutationEffect.getGeneAffected().equals(geneId)) {
				return true;
			}
		}
		return false;
	}

	public void addRawLine(String rawLine) {
		rawLines.add(rawLine);
	}

	public List<String> getRawLines() {
		return rawLines;
	}
	
	public String getFirstRawLine() {
		if (!rawLines.isEmpty()) {
			return rawLines.get(0);
		}
		return StringUtils.EMPTY;
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
