package org.coh.mattiyeh.datacruncher.model;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.StringUtils;
import org.coh.mattiyeh.datacruncher.genome.TriSeq;

public class Mutation implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 277334304446372552L;

	private String mutationId;
	private String donorId;
	private String specimenId;
	private String sampleId;
	private String matchedSampleId;
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
			String chr, int start, int end, String refBase, String mutBase, int totalReadCount,
			int mutantAlleleReadCount, String sequencingStrategy, String consequenceType, String geneAffected) {
		
		this(mutationId, donorId, specimenId, sampleId, matchedSampleId, chr, start, end, refBase, mutBase,
				totalReadCount, mutantAlleleReadCount, sequencingStrategy, consequenceType, geneAffected,
				new ArrayList<>());
	}

	public Mutation(String mutationId, String donorId, String specimenId, String sampleId, String matchedSampleId,
			String chr, int start, int end, String refBase, String mutBase, int totalReadCount,
			int mutantAlleleReadCount, String sequencingStrategy, String consequenceType, String geneAffected,
			List<String> rawLines) {
		super();
		this.mutationId = mutationId;
		this.donorId = donorId;
		this.specimenId = specimenId;
		this.sampleId = sampleId;
		this.matchedSampleId = matchedSampleId;
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

		this.rawLines = rawLines;
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

	public String getTriSeqWithMut() {
		if (triSeq == null) {
			return StringUtils.EMPTY;
		}
		String change = "T:" + triSeq.getPreBase() + "[" + triSeq.getRefBase() + ">" + mutBase + "]"
				+ triSeq.getPostBase();
		return change.toUpperCase();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((chr == null) ? 0 : chr.hashCode());
		result = prime * result + ((donorId == null) ? 0 : donorId.hashCode());
		result = prime * result + end;
		result = prime * result + (inCfsRegion ? 1231 : 1237);
		result = prime * result + (inPromoterRegion ? 1231 : 1237);
		result = prime * result + ((matchedSampleId == null) ? 0 : matchedSampleId.hashCode());
		result = prime * result + ((mutBase == null) ? 0 : mutBase.hashCode());
		result = prime * result + mutantAlleleReadCount;
		result = prime * result + ((mutationEffects == null) ? 0 : mutationEffects.hashCode());
		result = prime * result + ((mutationId == null) ? 0 : mutationId.hashCode());
		result = prime * result + ((rawLines == null) ? 0 : rawLines.hashCode());
		result = prime * result + ((refBase == null) ? 0 : refBase.hashCode());
		result = prime * result + ((sampleId == null) ? 0 : sampleId.hashCode());
		result = prime * result + ((sequencingStrategy == null) ? 0 : sequencingStrategy.hashCode());
		result = prime * result + ((specimenId == null) ? 0 : specimenId.hashCode());
		result = prime * result + start;
		result = prime * result + totalReadCount;
		result = prime * result + ((triSeq == null) ? 0 : triSeq.hashCode());
		return result;
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
		if (chr == null) {
			if (other.chr != null)
				return false;
		} else if (!chr.equals(other.chr))
			return false;
		if (donorId == null) {
			if (other.donorId != null)
				return false;
		} else if (!donorId.equals(other.donorId))
			return false;
		if (end != other.end)
			return false;
		if (inCfsRegion != other.inCfsRegion)
			return false;
		if (inPromoterRegion != other.inPromoterRegion)
			return false;
		if (matchedSampleId == null) {
			if (other.matchedSampleId != null)
				return false;
		} else if (!matchedSampleId.equals(other.matchedSampleId))
			return false;
		if (mutBase == null) {
			if (other.mutBase != null)
				return false;
		} else if (!mutBase.equals(other.mutBase))
			return false;
		if (mutantAlleleReadCount != other.mutantAlleleReadCount)
			return false;
		if (mutationEffects == null) {
			if (other.mutationEffects != null)
				return false;
		} else if (!mutationEffects.equals(other.mutationEffects))
			return false;
		if (mutationId == null) {
			if (other.mutationId != null)
				return false;
		} else if (!mutationId.equals(other.mutationId))
			return false;
		if (rawLines == null) {
			if (other.rawLines != null)
				return false;
		} else if (!rawLines.equals(other.rawLines))
			return false;
		if (refBase == null) {
			if (other.refBase != null)
				return false;
		} else if (!refBase.equals(other.refBase))
			return false;
		if (sampleId == null) {
			if (other.sampleId != null)
				return false;
		} else if (!sampleId.equals(other.sampleId))
			return false;
		if (sequencingStrategy == null) {
			if (other.sequencingStrategy != null)
				return false;
		} else if (!sequencingStrategy.equals(other.sequencingStrategy))
			return false;
		if (specimenId == null) {
			if (other.specimenId != null)
				return false;
		} else if (!specimenId.equals(other.specimenId))
			return false;
		if (start != other.start)
			return false;
		if (totalReadCount != other.totalReadCount)
			return false;
		if (triSeq == null) {
			if (other.triSeq != null)
				return false;
		} else if (!triSeq.equals(other.triSeq))
			return false;
		return true;
	}

}
