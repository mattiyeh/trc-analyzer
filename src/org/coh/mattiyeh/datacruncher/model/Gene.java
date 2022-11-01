package org.coh.mattiyeh.datacruncher.model;

public class Gene {

	private static final String PROTEIN_CODING = "protein_coding";
	
	String geneId;
	String geneSymbol;
	String chr;
	int start;
	int end;
	String geneType;

	public Gene(String geneId, String geneSymbol, String chr, int start, int end, String geneType) {
		super();
		this.geneId = geneId;
		this.geneSymbol = geneSymbol;
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.geneType = geneType;
	}

	public String getGeneId() {
		return geneId;
	}

	public String getGeneSymbol() {
		return geneSymbol;
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

	public String getGeneType() {
		return geneType;
	}

	public int getLength() {
		return end - start + 1;
	}

	public boolean isProteinCoding() {
		return PROTEIN_CODING.equals(geneType);
	}

	public boolean isNotProteinCoding() {
		return !isProteinCoding();
	}

}
