package org.coh.mattiyeh.datacruncher;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.coh.mattiyeh.datacruncher.model.Gene;
import org.coh.mattiyeh.datacruncher.model.GeneNaming;

public class GeneManager {

	private final File geneListFile = new File(Constants.WORKING_DIR, "data_files\\hg19_grch37_gene_list.txt");

	private Map<String, Gene> genesByGeneIds;
	private Map<String, Gene> genesByGeneSymbol;

	public GeneManager() throws NumberFormatException, IOException {
		// LOAD GENES
		genesByGeneIds = new HashMap<>();
		genesByGeneSymbol = new HashMap<>();

		Iterable<CSVRecord> records = CSVFormat.TDF.builder().setHeader().build().parse(new FileReader(geneListFile));
		for (CSVRecord tsvRecord : records) {
			String geneId = tsvRecord.get("gene_id");
			String geneSymbol = tsvRecord.get("gene_symbol");
			String chr = tsvRecord.get("chr");
			int start = Integer.parseInt(tsvRecord.get("start"));
			int end = Integer.parseInt(tsvRecord.get("end"));
			String geneType = tsvRecord.get("gene_type");

			Gene gene = new Gene(geneId, geneSymbol, chr, start, end, geneType);
			genesByGeneIds.put(geneId, gene);
			genesByGeneSymbol.put(geneSymbol, gene);
		}
	}

	public String getGeneSymbol(String geneId) {
		return genesByGeneIds.get(geneId).getGeneSymbol();
	}

	public Gene getGene(String name, GeneNaming gn) {
		switch (gn) {
		case ID:
			return genesByGeneIds.get(name);
		case SYMBOL:
			return genesByGeneSymbol.get(name);
		default:
			throw new IllegalArgumentException("Unexpected value: " + gn);
		}
	}

	public boolean geneInAssembly(String geneId) {
		return genesByGeneIds.containsKey(geneId);
	}

}
