package org.coh.mattiyeh.datacruncher;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.collections4.MapUtils;
import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.coh.mattiyeh.datacruncher.genome.GenomeExtractor;
import org.coh.mattiyeh.datacruncher.genome.TriSeq;
import org.coh.mattiyeh.datacruncher.math.Operator;
import org.coh.mattiyeh.datacruncher.model.DnaRange;
import org.coh.mattiyeh.datacruncher.model.Donor;
import org.coh.mattiyeh.datacruncher.model.Gene;
import org.coh.mattiyeh.datacruncher.model.GeneNaming;
import org.coh.mattiyeh.datacruncher.model.Mutation;
import org.coh.mattiyeh.datacruncher.model.MutationRange;
import org.coh.mattiyeh.datacruncher.model.MutationType;
import org.coh.mattiyeh.datacruncher.model.Sample;
import org.coh.mattiyeh.datacruncher.model.Specimen;

public class DataCruncher {

	private GeneManager gm;
	private GenomeExtractor ge;

	/**
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public void go() throws NumberFormatException, IOException {

		// Create output folder
		final String timestamp = new SimpleDateFormat("yyyy.MM.dd_HH.mm.ss").format(new Date());

		System.out.println("Initializing Gene Manager...");
		gm = new GeneManager();

		System.out.println("Initializing Genome Extractor...");
		ge = new GenomeExtractor(new File(Constants.DATA_FILES_DIR + "hg19.fa"));

		System.out.println();

		for (int i = 0; i < Constants.TUMOR_TYPES.length; i++) {

			String tumorType = Constants.TUMOR_TYPES[i];
			String tumorTypeDir = "icgc-dataset-" + tumorType;

			System.out.println("Starting " + tumorType);

			Path outputTimestampFolderPath = Paths.get(Constants.OUTPUT_DIR + "__" + timestamp, "DCO__" + timestamp + "_" + tumorType);
			Files.createDirectories(outputTimestampFolderPath);

			System.out.println("Reading donors...");
			Map<String, Donor> donors = readDonors(tumorTypeDir);

			System.out.println("Reading specimens...");
			readSpecimensAndSamples(tumorTypeDir, donors);

			System.out.println("Reading mutations...");
			readMutations(tumorTypeDir, donors);

			System.out.println("Reading expression data...");
			readExpressionLevels(tumorTypeDir, donors);

			System.out.println("Crunching numbers...");
			crunchNumbers(tumorType, tumorTypeDir, outputTimestampFolderPath, donors);

			System.out.println();
		}
		System.out.println("Done.");
	}

	private Map<String, Donor> readDonors(String tumorType) throws IOException {
		Map<String, Donor> donors = new HashMap<>();

		Path icgcDir = Paths.get(Constants.WORKING_ICGC_DIR);
		File donorsFile = icgcDir.resolve(tumorType).resolve(Constants.DONORS_FILENAME).toFile();

		Iterable<CSVRecord> records = CSVFormat.TDF.builder().setHeader().build()
				.parse(new InputStreamReader(new GZIPInputStream(new FileInputStream(donorsFile))));

		records.forEach(tsvRecord -> {
			String donorId = tsvRecord.get(Constants.ICGC_DONOR_ID);
			int survivalData = NumberUtils.toInt(tsvRecord.get("donor_survival_time"), -1);

			Donor donor = new Donor(donorId, survivalData);
			donors.put(donorId, donor);
		});

		return donors;
	}

	private void readSpecimensAndSamples(String tumorType, Map<String, Donor> donors) throws IOException {

		// Read in specimens
		Path icgcDir = Paths.get(Constants.WORKING_ICGC_DIR);
		File specimensFile = icgcDir.resolve(tumorType).resolve(Constants.SPECIMENS_FILENAME).toFile();

		Iterable<CSVRecord> records = CSVFormat.TDF.builder().setHeader().build()
				.parse(new InputStreamReader(new GZIPInputStream(new FileInputStream(specimensFile))));

		records.forEach(tsvRecord -> {
			String specimenId = tsvRecord.get(Constants.ICGC_SPECIMEN_ID);
			String donorId = tsvRecord.get(Constants.ICGC_DONOR_ID);
			
			String[] recordSpecimenType = tsvRecord.get(Constants.SPECIMEN_TYPE).split("-");
			String specimenType = recordSpecimenType[0].trim();
			String specimenSubType = recordSpecimenType[1].trim();
			
			String donorTreatmentType = tsvRecord.get(Constants.SPECIMEN_DONOR_TREATMENT_TYPE);

			// Reasons to skip specimens:
			// type is normal or xenograft or cell line
			// treatment is chemo or chemorads
			if (validSpecimenType(specimenType) && validSpecimenSubType(specimenSubType)
					&& validTreatmentType(donorTreatmentType)) {

				Specimen specimen = new Specimen(specimenId, donorId, specimenType, specimenSubType,
						donorTreatmentType);
				donors.get(donorId).addSpecimen(specimen);

			}
		});

		// Read in samples
		File samplesFile = icgcDir.resolve(tumorType).resolve(Constants.SAMPLES_FILENAME).toFile();
		records = CSVFormat.TDF.builder().setHeader().build()
				.parse(new InputStreamReader(new GZIPInputStream(new FileInputStream(samplesFile))));

		records.forEach(tsvRecord -> {
			String sampleId = tsvRecord.get(Constants.ICGC_SAMPLE_ID);
			String specimenId = tsvRecord.get(Constants.ICGC_SPECIMEN_ID);
			String donorId = tsvRecord.get(Constants.ICGC_DONOR_ID);
			Donor donor = donors.get(donorId);

			Sample sample = new Sample(sampleId, specimenId, donorId);

			// Need to check if specimen exists because it may not have passed QC from above
			if (donor.containsSpecimen(specimenId)) {
				donor.getSpecimen(specimenId).addSample(sample);
			}
		});

	}

	/**
	 * No specimens treated with chemo or rads
	 * 
	 * @param donorTreatmentType
	 * @return
	 */
	private boolean validTreatmentType(String donorTreatmentType) {
		return !(donorTreatmentType.equals(Constants.CHEMOTHERAPY)
				|| donorTreatmentType.equals(Constants.CHEMORADIATION));
	}

	/**
	 * No cell lines. No normals. No xenografts.
	 * 
	 * @param specimenType
	 * @return
	 */
	private boolean validSpecimenType(String specimenType) {
		return specimenType.equals(Constants.PRIMARY);
	}

	/**
	 * No blood derived
	 * 
	 * @param specimenSubType
	 * @return
	 */
	private boolean validSpecimenSubType(String specimenSubType) {
		return specimenSubType.equals(Constants.SOLID_TISSUE)
				|| specimenSubType.equals(Constants.ADDITIONAL_NEW_PRIMARY) || specimenSubType.equals(Constants.OTHER)
				|| specimenSubType.equals(Constants.LYMPH_NODE);
	}

	private void readMutations(String tumorType, Map<String, Donor> donors) throws IOException {

		// Read in promoter regions (data from UCSC)
		Map<String, List<DnaRange>> promoterRegionsByChr = readPromoterRegions();
		Map<String, List<DnaRange>> cfsRegionsByChr = readFragileSites();

		// Read in mutations
		Set<String> previousMutationKeys = new HashSet<>();

		CSVParser csvParser = openMutationFile(tumorType);

		for (CSVRecord tsvRecord : csvParser) {

			String mutationId = tsvRecord.get("icgc_mutation_id");
			String donorId = tsvRecord.get(Constants.ICGC_DONOR_ID);
			String specimenId = tsvRecord.get(Constants.ICGC_SPECIMEN_ID);
			String sampleId = tsvRecord.get(Constants.ICGC_SAMPLE_ID);
			String matchedSampleId = tsvRecord.get("matched_icgc_sample_id");
			String mutationType = tsvRecord.get("mutation_type");
			String chr = tsvRecord.get("chromosome");
			int start = Integer.parseInt(tsvRecord.get("chromosome_start"));
			int end = Integer.parseInt(tsvRecord.get("chromosome_end"));
			String refBase = tsvRecord.get("reference_genome_allele");
			String mutBase = tsvRecord.get("mutated_to_allele");
			int totalReadCount = NumberUtils.toInt(tsvRecord.get("total_read_count"));
			int mutantAlleleReadCount = NumberUtils.toInt(tsvRecord.get("mutant_allele_read_count"));
			String consequenceType = tsvRecord.get("consequence_type");
			String geneAffected = tsvRecord.get("gene_affected");
			String sequencingStrategy = tsvRecord.get(Constants.SEQUENCING_STRATEGY);

			Donor donor = donors.get(donorId);

			/*-
			 * Skip this mutation if:
			 * 
			 * 1. Gene affected column is blank (e.g. intergenic mutation)
			 * 2. Gene affected is not in GRch37 (hg19) Ensembl assembly
			 * 3. Gene is NOT protein coding
			 * 4. We threw out this specimen previously (see criteria in readSpecimensAndSamples)
			 * 5. Not WGS
			 */
			if (StringUtils.isBlank(geneAffected) || !gm.geneInAssembly(geneAffected)
					|| gm.getGene(geneAffected, GeneNaming.ID).isNotProteinCoding()
					|| !donor.containsSpecimen(specimenId) || !Constants.WGS.equals(sequencingStrategy)) {
				continue;
			}

			// Can't just use mutation ID from the data because it is NOT unique (e.g. same
			// mutation id used in different samples)
			String mutationKey = new StringBuffer().append(mutationId).append(donorId).append(specimenId)
					.append(sampleId).append(chr).append(start).append(end).append(refBase).append(mutBase).toString();
			
			Sample sample = donor.getSpecimen(specimenId).getSample(sampleId);

			// Check to see if we've seen this mutation for this sample before (ie.
			// different consequence or gene affected)
			if (previousMutationKeys.contains(mutationKey)) {

				Mutation oldMutation = sample.getMutation(mutationId);

				// Check to see if the previously seen mutation already covered this gene
				if (!oldMutation.affectsGene(geneAffected)) {
					oldMutation.addMutationEffect(consequenceType, geneAffected);
					oldMutation.addRawLine(StringUtils.join(tsvRecord.toList(), '\t'));
				}

			} else {

				Mutation newMutation = new Mutation(mutationId, donorId, specimenId, sampleId, matchedSampleId,
						mutationType, chr, start, end, refBase, mutBase, totalReadCount, mutantAlleleReadCount,
						sequencingStrategy, consequenceType, geneAffected);

				newMutation.addRawLine(StringUtils.join(tsvRecord.toList(), '\t'));

				if (Constants.SBS.equals(mutationType)) {
					TriSeq triSeq = ge.getTrinucleotideContext(chr, start);
					newMutation.setTriSeq(triSeq);
				}

				// Retrieve promoter ranges in null-safe manner in case of unusual chr (e.g.
				// "MT")
				List<DnaRange> promoterRanges = MapUtils.getObject(promoterRegionsByChr, chr, new ArrayList<>());
				for (DnaRange promoterRange : promoterRanges) {
					if (promoterRange.overlaps(start, end)) {
						newMutation.setInPromoterRegion(true);

						// No need to keep checking the rest of the promoter ranges
						break;
					}
				}

				// Retrieve CFS ranges in null-safe manner in case of unusual chr (e.g. "MT")
				List<DnaRange> cfsRanges = MapUtils.getObject(cfsRegionsByChr, chr, new ArrayList<>());
				for (DnaRange cfsRange : cfsRanges) {
					if (cfsRange.overlaps(start, end)) {
						newMutation.setInCfsRegion(true);

						// No need to keep checking the rest of the CFS ranges
						break;
					}
				}

				previousMutationKeys.add(mutationKey);

				sample.addMutation(newMutation);

			}
		}
		csvParser.close();
	}

	/**
	 * @param specimens
	 * @throws IOException
	 */
	private void readExpressionLevels(String tumorType, Map<String, Donor> donors) throws IOException {

		Path icgcDir = Paths.get(Constants.WORKING_ICGC_DIR);
		File expressionFile = icgcDir.resolve(tumorType).resolve("exp_seq.tsv.gz").toFile();

		// check to make sure we have expression data for this tissue type
		if (!expressionFile.exists()) {
			System.out.println("WARNING: No expression data found");
			return;
		}

		Iterable<CSVRecord> records = CSVFormat.TDF.builder().setHeader().build()
				.parse(new InputStreamReader(new GZIPInputStream(new FileInputStream(expressionFile))));

		for (CSVRecord tsvRecord : records) {
			String donorId = tsvRecord.get(Constants.ICGC_DONOR_ID);
			String specimenId = tsvRecord.get(Constants.ICGC_SPECIMEN_ID);
			String sampleId = tsvRecord.get(Constants.ICGC_SAMPLE_ID);
			String geneId = tsvRecord.get("gene_id");
			float normalizedReadCount = NumberUtils.toFloat(tsvRecord.get("normalized_read_count"));

			// Check to make sure that geneId is the actual ENSG id and not the gene symbol
			if (!geneId.startsWith("ENSG")) {
				// Find the corresponding gene id
				Gene gene = gm.getGene(geneId, GeneNaming.SYMBOL);
				if (gene == null) {
					// This gene is not in our assembly. Skip.
					continue;
				}
				geneId = gene.getGeneId();
			}

			if (donors.get(donorId).containsSpecimen(specimenId)) {
				Sample sample = donors.get(donorId).getSpecimen(specimenId).getSample(sampleId);
				sample.addGeneNormExpressionLevel(geneId, normalizedReadCount);
			}
		}

	}

	private void crunchNumbers(String tumorType, String tumorTypeDir, Path outputFolderPath, Map<String, Donor> donors)
			throws IOException {

		final int highCutoff = 75;
		final int lowCutoff = 25;
		final int zeroCutoff = 0;
		
		Path promoterPath = outputFolderPath.resolve("promoter");
		Path promoterSbsPath = promoterPath.resolve("sbs");
		Path promoterIndelPath = promoterPath.resolve("indel");
		Path promoterMbsPath = promoterPath.resolve("mbs");
		Files.createDirectories(promoterSbsPath);
		Files.createDirectories(promoterIndelPath);
		Files.createDirectories(promoterMbsPath);
		
		Path cfsPath = outputFolderPath.resolve("cfs");
		Path cfsSbsPath = cfsPath.resolve("sbs");
		Path cfsIndelPath = cfsPath.resolve("indel");
		Path cfsMbsPath = cfsPath.resolve("mbs");
		Files.createDirectories(cfsSbsPath);
		Files.createDirectories(cfsIndelPath);
		Files.createDirectories(cfsMbsPath);
		
		Path metadataPath = outputFolderPath.resolve(tumorType + "_metadata.tsv");
		
		Path mutationPath = outputFolderPath.resolve(tumorType + "_mutations.tsv");
		
		// PROMOTER
		
		Path nonPromoterMutationPath = promoterPath.resolve(tumorType + "_non_promoter_mutations.tsv");
		Path promoterMutationPath = promoterPath.resolve(tumorType + "_promoter_mutations.tsv");
		
		// PROMOTER-SBS
		
		Path nonPromoterSbsMutationPath = promoterSbsPath.resolve(tumorType + "_non_promoter_sbs_mutations.tsv");
		Path promoterSbsMutationPath = promoterSbsPath.resolve(tumorType + "_promoter_sbs_mutations.tsv");		
		Path promoterSbsMutationHighPath = promoterSbsPath.resolve(tumorType + "_promoter_sbs_mutations_high.tsv");		
		Path promoterSbsMutationMidPath = promoterSbsPath.resolve(tumorType + "_promoter_sbs_mutations_mid.tsv");		
		Path promoterSbsMutationLowPath = promoterSbsPath.resolve(tumorType + "_promoter_sbs_mutations_low.tsv");		
		Path promoterSbsMutationZeroPath = promoterSbsPath.resolve(tumorType + "_promoter_sbs_mutations_zero.tsv");		
		
		// PROMOTER-INDEL
		
		Path nonPromoterIndelMutationPath = promoterIndelPath.resolve(tumorType + "_non_promoter_indel_mutations.tsv");
		Path promoterIndelMutationPath = promoterIndelPath.resolve(tumorType + "_promoter_indel_mutations.tsv");		
		Path promoterIndelMutationHighPath = promoterIndelPath.resolve(tumorType + "_promoter_indel_mutations_high.tsv");		
		Path promoterIndelMutationMidPath = promoterIndelPath.resolve(tumorType + "_promoter_indel_mutations_mid.tsv");		
		Path promoterIndelMutationLowPath = promoterIndelPath.resolve(tumorType + "_promoter_indel_mutations_low.tsv");		
		Path promoterIndelMutationZeroPath = promoterIndelPath.resolve(tumorType + "_promoter_indel_mutations_zero.tsv");		
			
		// PROMOTER-MBS
		
		Path nonPromoterMbsMutationPath = promoterMbsPath.resolve(tumorType + "_non_promoter_mbs_mutations.tsv");
		Path promoterMbsMutationPath = promoterMbsPath.resolve(tumorType + "_promoter_mbs_mutations.tsv");		
		Path promoterMbsMutationHighPath = promoterMbsPath.resolve(tumorType + "_promoter_mbs_mutations_high.tsv");		
		Path promoterMbsMutationMidPath = promoterMbsPath.resolve(tumorType + "_promoter_mbs_mutations_mid.tsv");		
		Path promoterMbsMutationLowPath = promoterMbsPath.resolve(tumorType + "_promoter_mbs_mutations_low.tsv");		
		Path promoterMbsMutationZeroPath = promoterMbsPath.resolve(tumorType + "_promoter_mbs_mutations_zero.tsv");		
		
		// CFS
		
		Path nonCfsMutationsPath = cfsPath.resolve(tumorType + "_non_cfs_mutations.tsv");
		Path cfsMutationsPath = cfsPath.resolve(tumorType + "_cfs_mutations.tsv");
		
		// CFS-SBS
		
		Path nonCfsSbsMutationsPath = cfsSbsPath.resolve(tumorType + "_non_cfs_sbs_mutations.tsv");
		Path cfsSbsMutationsPath = cfsSbsPath.resolve(tumorType + "_cfs_sbs_mutations.tsv");
		
		// CFS-INDEL
		
		Path nonCfsIndelMutationsPath = cfsIndelPath.resolve(tumorType + "_non_cfs_indel_mutations.tsv");
		Path cfsIndelMutationsPath = cfsIndelPath.resolve(tumorType + "_cfs_indel_mutations.tsv");
		
		// CFS-MBS
		
		Path nonCfsMbsMutationsPath = cfsMbsPath.resolve(tumorType + "_non_cfs_mbs_mutations.tsv");
		Path cfsMbsMutationsPath = cfsMbsPath.resolve(tumorType + "_cfs_mbs_mutations.tsv");
		
		try (BufferedWriter metadataBw = Files.newBufferedWriter(metadataPath);
				
				BufferedWriter mutsBw = Files.newBufferedWriter(mutationPath);
				
				BufferedWriter nonPromMutsBw = Files.newBufferedWriter(nonPromoterMutationPath);
				BufferedWriter promMutsBw = Files.newBufferedWriter(promoterMutationPath);
				
				BufferedWriter nonPromSbsMutsBw = Files.newBufferedWriter(nonPromoterSbsMutationPath);
				BufferedWriter promSbsMutsBw = Files.newBufferedWriter(promoterSbsMutationPath);
				BufferedWriter promSbsMutsHighBw = Files.newBufferedWriter(promoterSbsMutationHighPath);
				BufferedWriter promSbsMutsMidBw = Files.newBufferedWriter(promoterSbsMutationMidPath);
				BufferedWriter promSbsMutsLowBw = Files.newBufferedWriter(promoterSbsMutationLowPath);
				BufferedWriter promSbsMutsZeroBw = Files.newBufferedWriter(promoterSbsMutationZeroPath);
				
				BufferedWriter nonPromIndelMutsBw = Files.newBufferedWriter(nonPromoterIndelMutationPath);
				BufferedWriter promIndelMutsBw = Files.newBufferedWriter(promoterIndelMutationPath);
				BufferedWriter promIndelMutsHighBw = Files.newBufferedWriter(promoterIndelMutationHighPath);
				BufferedWriter promIndelMutsMidBw = Files.newBufferedWriter(promoterIndelMutationMidPath);
				BufferedWriter promIndelMutsLowBw = Files.newBufferedWriter(promoterIndelMutationLowPath);
				BufferedWriter promIndelMutsZeroBw = Files.newBufferedWriter(promoterIndelMutationZeroPath);
				
				BufferedWriter nonPromMbsMutsBw = Files.newBufferedWriter(nonPromoterMbsMutationPath);
				BufferedWriter promMbsMutsBw = Files.newBufferedWriter(promoterMbsMutationPath);
				BufferedWriter promMbsMutsHighBw = Files.newBufferedWriter(promoterMbsMutationHighPath);
				BufferedWriter promMbsMutsMidBw = Files.newBufferedWriter(promoterMbsMutationMidPath);
				BufferedWriter promMbsMutsLowBw = Files.newBufferedWriter(promoterMbsMutationLowPath);
				BufferedWriter promMbsMutsZeroBw = Files.newBufferedWriter(promoterMbsMutationZeroPath);
				
				BufferedWriter nonCfsMutsBw = Files.newBufferedWriter(nonCfsMutationsPath);
				BufferedWriter cfsMutsBw = Files.newBufferedWriter(cfsMutationsPath);
				
				BufferedWriter nonCfsSbsMutsBw = Files.newBufferedWriter(nonCfsSbsMutationsPath);
				BufferedWriter cfsSbsMutsBw = Files.newBufferedWriter(cfsSbsMutationsPath);
				
				BufferedWriter nonCfsIndelMutsBw = Files.newBufferedWriter(nonCfsIndelMutationsPath);
				BufferedWriter cfsIndelMutsBw = Files.newBufferedWriter(cfsIndelMutationsPath);
				
				BufferedWriter nonCfsMbsMutsBw = Files.newBufferedWriter(nonCfsMbsMutationsPath);
				BufferedWriter cfsMbsMutsBw = Files.newBufferedWriter(cfsMbsMutationsPath);) {

			List<String> headerItems = new ArrayList<>();
			headerItems.add("donor_id");
			headerItems.add("specimens");
			headerItems.add("samples");
			headerItems.add("spec_w_muts");
			headerItems.add("spec_w_exp");
			headerItems.add("spec_w_both");
			
			headerItems.add("all_muts");
			
			headerItems.add("non_prom_muts");
			headerItems.add("prom_muts");
			
			headerItems.add("non_prom_sbs_muts");
			headerItems.add("prom_sbs_muts");
			headerItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff);
			headerItems.add("prom_sbs_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
			headerItems.add("prom_sbs_muts_w_exp_lte_" + lowCutoff);
			headerItems.add("prom_sbs_muts_w_zero");
			
			headerItems.add("non_prom_indel_muts");
			headerItems.add("prom_indel_muts");
			headerItems.add("prom_indel_muts_w_exp_gte_" + highCutoff);
			headerItems.add("prom_indel_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
			headerItems.add("prom_indel_muts_w_exp_lte_" + lowCutoff);
			headerItems.add("prom_indel_muts_w_zero");
			
			headerItems.add("non_prom_mbs_muts");
			headerItems.add("prom_mbs_muts");
			headerItems.add("prom_mbs_muts_w_exp_gte_" + highCutoff);
			headerItems.add("prom_mbs_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
			headerItems.add("prom_mbs_muts_w_exp_lte_" + lowCutoff);
			headerItems.add("prom_mbs_muts_w_zero");
			
			headerItems.add("non_cfs_muts");
			headerItems.add("cfs_muts");
			
			headerItems.add("non_cfs_sbs_muts");
			headerItems.add("cfs_sbs_muts");
			
			headerItems.add("non_cfs_indel_muts");
			headerItems.add("cfs_indel_muts");
			
			headerItems.add("non_cfs_mbs_muts");
			headerItems.add("cfs_mbs_muts");
						
			metadataBw.write(StringUtils.join(headerItems, '\t'));
			metadataBw.newLine();

			// Grab mutation file header
			CSVParser csvParser = openMutationFile(tumorTypeDir);
			List<String> headerElements = new ArrayList<>(csvParser.getHeaderNames());
			csvParser.close();
			headerElements.add("Ref_Tri");
			headerElements.add("genes_affected");
			headerElements.add("all_genes_affected");
			String header = StringUtils.join(headerElements, '\t');
			
			writeHeader(mutsBw, header);
			
			writeHeader(nonPromMutsBw, header);
			writeHeader(promMutsBw, header);
			
			writeHeader(nonPromSbsMutsBw, header);
			writeHeader(promSbsMutsBw, header);
			writeHeader(promSbsMutsHighBw, header);
			writeHeader(promSbsMutsMidBw, header);
			writeHeader(promSbsMutsLowBw, header);
			writeHeader(promSbsMutsZeroBw, header);
			
			writeHeader(nonPromIndelMutsBw, header);
			writeHeader(promIndelMutsBw, header);
			writeHeader(promIndelMutsHighBw, header);
			writeHeader(promIndelMutsMidBw, header);
			writeHeader(promIndelMutsLowBw, header);
			writeHeader(promIndelMutsZeroBw, header);
			
			writeHeader(nonPromMbsMutsBw, header);
			writeHeader(promMbsMutsBw, header);
			writeHeader(promMbsMutsHighBw, header);
			writeHeader(promMbsMutsMidBw, header);
			writeHeader(promMbsMutsLowBw, header);
			writeHeader(promMbsMutsZeroBw, header);
			
			writeHeader(nonCfsMutsBw, header);
			writeHeader(cfsMutsBw, header);
			
			writeHeader(nonCfsSbsMutsBw, header);
			writeHeader(cfsSbsMutsBw, header);
			
			writeHeader(nonCfsIndelMutsBw, header);
			writeHeader(cfsIndelMutsBw, header);
			
			writeHeader(nonCfsMbsMutsBw, header);
			writeHeader(cfsMbsMutsBw, header);

			for (Map.Entry<String, Donor> entry : donors.entrySet()) {
				String donorId = entry.getKey();
				Donor donor = entry.getValue();

				// Only work with donors that have a specimen that has BOTH (based on discussion
				// with Raoof on 5/8/2022)
				if (!donor.hasMutationAndExpressionData()) {
					continue;
				}

				Set<Mutation> mutations = donor.getMutations(MutationRange.NONE, MutationType.ALL);
				
				// PROMOTER
				
				Set<Mutation> nonPromoterMutations = donor.getMutations(MutationRange.NONPROMOTER, MutationType.ALL);
				Set<Mutation> promoterMutations = donor.getMutations(MutationRange.PROMOTER, MutationType.ALL);
				
				// PROMOTER-SBS
				
				Set<Mutation> nonPromoterSbsMutations = donor.getMutations(MutationRange.NONPROMOTER, MutationType.SBS);
				Set<Mutation> promoterSbsMutations = donor.getMutations(MutationRange.PROMOTER, MutationType.SBS);
				// Get promoter SBS mutations in high expressed genes (highCutoff <= HERE)
				Set<Mutation> promoterSbsMutationsInHighExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.SBS, Operator.GREATERTHANOREQUAL, highCutoff);
				// Get promoter SBS mutations in mid expressed genes (lowCutoff < HERE < highCutoff)
				Set<Mutation> promoterSbsMutationsInMidExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.SBS, Operator.GREATERTHAN, lowCutoff, Operator.LESSTHAN, highCutoff);
				// Get promoter SBS mutations in low expressed genes (HERE <= lowCutoff)
				Set<Mutation> promoterSbsMutationsInLowExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.SBS, Operator.LESSTHANOREQUAL, lowCutoff);
				// Get promoter SBS mutations in NON-expressed genes (HERE == zeroCutoff)
				Set<Mutation> promoterSbsMutationsInZeroExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.SBS, Operator.EQUALTO, zeroCutoff);
				
				// PROMOTER-INDEL
				
				Set<Mutation> nonPromoterIndelMutations = donor.getIndelMutations(MutationRange.NONPROMOTER);
				Set<Mutation> promoterIndelMutations = donor.getIndelMutations(MutationRange.PROMOTER);
				// Get promoter INDEL mutations in high expressed genes (highCutoff <= HERE)
				Set<Mutation> promoterIndelMutationsInHighExpressedGenes = donor.getIndelMutationsInExpressedGenes(MutationRange.PROMOTER, Operator.GREATERTHANOREQUAL, highCutoff);
				// Get promoter INDEL mutations in mid expressed genes (lowCutoff < HERE < highCutoff)
				Set<Mutation> promoterIndelMutationsInMidExpressedGenes = donor.getIndelMutationsInExpressedGenes(MutationRange.PROMOTER, Operator.GREATERTHAN, lowCutoff, Operator.LESSTHAN, highCutoff);
				// Get promoter INDEL mutations in low expressed genes (HERE <= lowCutoff)
				Set<Mutation> promoterIndelMutationsInLowExpressedGenes = donor.getIndelMutationsInExpressedGenes(MutationRange.PROMOTER, Operator.LESSTHANOREQUAL, lowCutoff);
				// Get promoter INDEL mutations in NON-expressed genes (HERE == zeroCutoff)
				Set<Mutation> promoterIndelMutationsInZeroExpressedGenes = donor.getIndelMutationsInExpressedGenes(MutationRange.PROMOTER, Operator.EQUALTO, zeroCutoff);

				// PROMOTER-MBS
				
				Set<Mutation> nonPromoterMbsMutations = donor.getMutations(MutationRange.NONPROMOTER, MutationType.MBS);
				Set<Mutation> promoterMbsMutations = donor.getMutations(MutationRange.PROMOTER, MutationType.MBS);
				// Get promoter SBS mutations in high expressed genes (highCutoff <= HERE)
				Set<Mutation> promoterMbsMutationsInHighExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.MBS, Operator.GREATERTHANOREQUAL, highCutoff);
				// Get promoter SBS mutations in mid expressed genes (lowCutoff < HERE < highCutoff)
				Set<Mutation> promoterMbsMutationsInMidExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.MBS, Operator.GREATERTHAN, lowCutoff, Operator.LESSTHAN, highCutoff);
				// Get promoter SBS mutations in low expressed genes (HERE <= lowCutoff)
				Set<Mutation> promoterMbsMutationsInLowExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.MBS, Operator.LESSTHANOREQUAL, lowCutoff);
				// Get promoter SBS mutations in NON-expressed genes (HERE == zeroCutoff)
				Set<Mutation> promoterMbsMutationsInZeroExpressedGenes = donor.getMutationsInExpressedGenes(MutationRange.PROMOTER, MutationType.MBS, Operator.EQUALTO, zeroCutoff);
				
				// CFS
				
				Set<Mutation> nonCfsMutations = donor.getMutations(MutationRange.NONCFS, MutationType.ALL);
				Set<Mutation> cfsMutations = donor.getMutations(MutationRange.CFS, MutationType.ALL);
				
				// CFS-SBS
				
				Set<Mutation> nonCfsSbsMutations = donor.getMutations(MutationRange.NONCFS, MutationType.SBS);
				Set<Mutation> cfsSbsMutations = donor.getMutations(MutationRange.CFS, MutationType.SBS);

				// CFS-INDEL
				
				Set<Mutation> nonCfsIndelMutations = donor.getIndelMutations(MutationRange.NONCFS);
				Set<Mutation> cfsIndelMutations = donor.getIndelMutations(MutationRange.CFS);
				
				// CFS-MBS
				
				Set<Mutation> nonCfsMbsMutations = donor.getMutations(MutationRange.NONCFS, MutationType.MBS);
				Set<Mutation> cfsMbsMutations = donor.getMutations(MutationRange.CFS, MutationType.MBS);
				
				metadataBw.write(donorId);
				metadataBw.write("\t" + donor.getNumSpecimens());
				metadataBw.write("\t" + donor.getNumSamples());
				metadataBw.write("\t" + donor.getNumSpecimensWithMutationData());
				metadataBw.write("\t" + donor.getNumSpecimensWithExpressionData());
				metadataBw.write("\t" + donor.getNumSpecimensWithBoth());
				
				metadataBw.write("\t" + mutations.size());
				
				metadataBw.write("\t" + nonPromoterMutations.size());
				metadataBw.write("\t" + promoterMutations.size());
				
				metadataBw.write("\t" + nonPromoterSbsMutations.size());
				metadataBw.write("\t" + promoterSbsMutations.size());
				metadataBw.write("\t" + promoterSbsMutationsInHighExpressedGenes.size());
				metadataBw.write("\t" + promoterSbsMutationsInMidExpressedGenes.size());
				metadataBw.write("\t" + promoterSbsMutationsInLowExpressedGenes.size());
				metadataBw.write("\t" + promoterSbsMutationsInZeroExpressedGenes.size());
				
				metadataBw.write("\t" + nonPromoterIndelMutations.size());
				metadataBw.write("\t" + promoterIndelMutations.size());
				metadataBw.write("\t" + promoterIndelMutationsInHighExpressedGenes.size());
				metadataBw.write("\t" + promoterIndelMutationsInMidExpressedGenes.size());
				metadataBw.write("\t" + promoterIndelMutationsInLowExpressedGenes.size());
				metadataBw.write("\t" + promoterIndelMutationsInZeroExpressedGenes.size());
				
				metadataBw.write("\t" + nonPromoterMbsMutations.size());
				metadataBw.write("\t" + promoterMbsMutations.size());
				metadataBw.write("\t" + promoterMbsMutationsInHighExpressedGenes.size());
				metadataBw.write("\t" + promoterMbsMutationsInMidExpressedGenes.size());
				metadataBw.write("\t" + promoterMbsMutationsInLowExpressedGenes.size());
				metadataBw.write("\t" + promoterMbsMutationsInZeroExpressedGenes.size());
				
				metadataBw.write("\t" + nonCfsMutations.size());
				metadataBw.write("\t" + cfsMutations.size());
				
				metadataBw.write("\t" + nonCfsSbsMutations.size());
				metadataBw.write("\t" + cfsSbsMutations.size());
				
				metadataBw.write("\t" + nonCfsIndelMutations.size());
				metadataBw.write("\t" + cfsIndelMutations.size());
				
				metadataBw.write("\t" + nonCfsMbsMutations.size());
				metadataBw.write("\t" + cfsMbsMutations.size());
				
				metadataBw.newLine();
				
				writeMutationsToFiles(mutations, mutsBw);
				
				// PROMOTER
				
				writeMutationsToFiles(nonPromoterMutations, nonPromMutsBw);
				writeMutationsToFiles(promoterMutations, promMutsBw);

				// PROMOTER-SBS
								
				writeMutationsToFiles(nonPromoterSbsMutations, nonPromSbsMutsBw);
				writeMutationsToFiles(promoterSbsMutations, promSbsMutsBw);
				writeMutationsToFiles(promoterSbsMutationsInHighExpressedGenes, promSbsMutsHighBw);
				writeMutationsToFiles(promoterSbsMutationsInMidExpressedGenes, promSbsMutsMidBw);
				writeMutationsToFiles(promoterSbsMutationsInLowExpressedGenes, promSbsMutsLowBw);
				writeMutationsToFiles(promoterSbsMutationsInZeroExpressedGenes, promSbsMutsZeroBw);
				
				// PROMOTER-INDEL
				
				writeMutationsToFiles(nonPromoterIndelMutations, nonPromIndelMutsBw);
				writeMutationsToFiles(promoterIndelMutations, promIndelMutsBw);
				writeMutationsToFiles(promoterIndelMutationsInHighExpressedGenes, promIndelMutsHighBw);
				writeMutationsToFiles(promoterIndelMutationsInMidExpressedGenes, promIndelMutsMidBw);
				writeMutationsToFiles(promoterIndelMutationsInLowExpressedGenes, promIndelMutsLowBw);
				writeMutationsToFiles(promoterIndelMutationsInZeroExpressedGenes, promIndelMutsZeroBw);
				
				// PROMOTER-MBS
				
				writeMutationsToFiles(nonPromoterMbsMutations, nonPromMbsMutsBw);
				writeMutationsToFiles(promoterMbsMutations, promMbsMutsBw);
				writeMutationsToFiles(promoterMbsMutationsInHighExpressedGenes, promMbsMutsHighBw);
				writeMutationsToFiles(promoterMbsMutationsInMidExpressedGenes, promMbsMutsMidBw);
				writeMutationsToFiles(promoterMbsMutationsInLowExpressedGenes, promMbsMutsLowBw);
				writeMutationsToFiles(promoterMbsMutationsInZeroExpressedGenes, promMbsMutsZeroBw);
								
				// CFS

				writeMutationsToFiles(nonCfsMutations, nonCfsMutsBw);
				writeMutationsToFiles(cfsMutations, cfsMutsBw);
				
				// CFS-SBS

				writeMutationsToFiles(nonCfsSbsMutations, nonCfsSbsMutsBw);
				writeMutationsToFiles(cfsSbsMutations, cfsSbsMutsBw);

				// CFS-INDEL

				writeMutationsToFiles(nonCfsIndelMutations, nonCfsIndelMutsBw);
				writeMutationsToFiles(cfsIndelMutations, cfsIndelMutsBw);

				// CFS-MBS

				writeMutationsToFiles(nonCfsMbsMutations, nonCfsMbsMutsBw);
				writeMutationsToFiles(cfsMbsMutations, cfsMbsMutsBw);

			}
		}
	}

	private void writeHeader(BufferedWriter bw, String header) throws IOException {
		bw.write(header);
		bw.newLine();
	}

	private void writeMutationsToFiles(Set<Mutation> promoterMutations, BufferedWriter bw) throws IOException {
		for (Mutation mutation : promoterMutations) {
			
			// Write the main info only from the first version of the mutation. No need for
			// duplicate mutation lines with different affected genes
			bw.write(mutation.getFirstRawLine());
			
			// Add three-base data to end of line
			bw.write("\t" + mutation.getTriSeqWithMut());
			
			// Concatenate genes affected and put it in last column
			bw.write("\t" + mutation.getNumGeneIdsAffected());
			bw.write("\t" + StringUtils.join(mutation.getGeneIdsAffected(), ';'));
			
			bw.newLine();
		}
	}

	private CSVParser openMutationFile(String tumorTypeDir) throws IOException {
		Path icgcDir = Paths.get(Constants.WORKING_ICGC_DIR);
		File mutationsFile = icgcDir.resolve(tumorTypeDir).resolve(Constants.MUTATIONS_FILENAME).toFile();

		CSVFormat csvFormat = CSVFormat.TDF.builder().setHeader().build();
		return CSVParser.parse(new InputStreamReader(new GZIPInputStream(new FileInputStream(mutationsFile))),
				csvFormat);
	}

	private Map<String, List<DnaRange>> readPromoterRegions() throws IOException {
		// Read in all promoter regions (from UCSC upstream1000.fa.gz file)
		// http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/
		/*
		 * Sequences 1000 bases upstream of annotated transcription starts of RefSeq
		 * genes with annotated 5' UTRs. This file is updated weekly so it might be
		 * slightly out of sync with the RefSeq data which is updated daily for most
		 * assemblies.
		 */
		File promoterRegionsFile = new File(Constants.DATA_FILES_DIR + "promoter_regions.txt");
		final String promoterRegex = "chr(X|Y|\\d+):(\\d+)-(\\d+)";
		final Pattern promoterPattern = Pattern.compile(promoterRegex, Pattern.MULTILINE);
		Map<String, List<DnaRange>> promoterRegionsByChr = new HashMap<>();
		String st;
		try (BufferedReader br = new BufferedReader(new FileReader(promoterRegionsFile))) {
			while ((st = br.readLine()) != null) {
				final Matcher matcher = promoterPattern.matcher(st);
				while (matcher.find()) {

					String chr = matcher.group(1);
					int promoterRangeStart = Integer.parseInt(matcher.group(2));
					int promoterRangeEnd = Integer.parseInt(matcher.group(3));

					promoterRegionsByChr.putIfAbsent(chr, new ArrayList<>());

					List<DnaRange> promoterRegions = promoterRegionsByChr.get(chr);
					DnaRange promoterRegion = new DnaRange(chr, promoterRangeStart, promoterRangeEnd);
					promoterRegions.add(promoterRegion);

				}
			}
		}
		return promoterRegionsByChr;
	}

	private Map<String, List<DnaRange>> readFragileSites() throws IOException {
		// Read in all fragile sites (from Table 1; Ma et al. 2012 Int J Mol Sci)
		// Ma K, Qiu L, Mrasek K, Zhang J, Liehr T, Quintana LG, Li Z. Common fragile
		// sites: genomic hotspots of DNA damage and carcinogenesis. Int J Mol Sci.
		// 2012;13(9):11974-11999. doi: 10.3390/ijms130911974. Epub 2012 Sep 20. PMID:
		// 23109895; PMCID: PMC3472787.
		
		// Cytoband locations from UCSC cytoband.txt.gz hg19 merged into cfs_sites.tsv

		Map<String, List<DnaRange>> fragileSitesByChr = new HashMap<>();

		File fragileSitesFile = new File(Constants.DATA_FILES_DIR + "cfs_sites.tsv");
		Iterable<CSVRecord> records = CSVFormat.TDF.builder().setHeader().build()
				.parse(new FileReader(fragileSitesFile));

		records.forEach(tsvRecord -> {
			String chr = tsvRecord.get("chr");
			int start = Integer.parseInt(tsvRecord.get("start"));
			int end = Integer.parseInt(tsvRecord.get("end"));

			fragileSitesByChr.putIfAbsent(chr, new ArrayList<>());

			DnaRange cfs = new DnaRange(chr, start, end);
			List<DnaRange> fragileSites = fragileSitesByChr.get(chr);
			fragileSites.add(cfs);
		});

		return fragileSitesByChr;
	}

}
