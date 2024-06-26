package org.coh.mattiyeh.datacruncher.model;

import java.io.Serializable;

public class MutationEffect implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7408995583952012439L;

	String consequenceType;
	String geneAffected;
	String transcriptAffected;
	boolean proteinCoding;
	String rawLine;

	public MutationEffect(String consequenceType, String geneAffected, String transcriptAffected, boolean proteinCoding, String rawLine) {
		super();
		this.consequenceType = consequenceType;
		this.geneAffected = geneAffected;
		this.transcriptAffected = transcriptAffected;
		this.proteinCoding = proteinCoding;
		this.rawLine = rawLine;
	}

	public String getConsequenceType() {
		return consequenceType;
	}

	public String getGeneAffected() {
		return geneAffected;
	}
	
	public String getTranscriptAffected() {
		return transcriptAffected;
	}
	
	public boolean isProteinCoding() {
		return proteinCoding;
	}
	
	public String getRawLine() {
		return rawLine;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((consequenceType == null) ? 0 : consequenceType.hashCode());
		result = prime * result + ((geneAffected == null) ? 0 : geneAffected.hashCode());
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
		MutationEffect other = (MutationEffect) obj;
		if (consequenceType == null) {
			if (other.consequenceType != null)
				return false;
		} else if (!consequenceType.equals(other.consequenceType))
			return false;
		if (geneAffected == null) {
			if (other.geneAffected != null)
				return false;
		} else if (!geneAffected.equals(other.geneAffected))
			return false;
		return true;
	}

}
