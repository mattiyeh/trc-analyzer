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

			Path outputTimestampFolderPath = Paths.get(Constants.WORKING_DIR,
					"DataCruncherOutput__" + timestamp + "_" + tumorType);
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
		File donorsFile = icgcDir.resolve(tumorType).resolve("donor.tsv.gz").toFile();

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
		File specimensFile = icgcDir.resolve(tumorType).resolve("specimen.tsv.gz").toFile();

		Iterable<CSVRecord> records = CSVFormat.TDF.builder().setHeader().build()
				.parse(new InputStreamReader(new GZIPInputStream(new FileInputStream(specimensFile))));

		records.forEach(tsvRecord -> {
			String specimenId = tsvRecord.get(Constants.ICGC_SPECIMEN_ID);
			String donorId = tsvRecord.get(Constants.ICGC_DONOR_ID);
			String specimenType = tsvRecord.get("specimen_type").split("-")[0].trim();
			String specimenSubType = tsvRecord.get("specimen_type").split("-")[1].trim();
			String donorTreatmentType = tsvRecord.get("specimen_donor_treatment_type");

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
		File samplesFile = icgcDir.resolve(tumorType).resolve("sample.tsv.gz").toFile();
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

			// Can't use mutation ID from the data because it is NOT unique (e.g. same
			// mutation id used in different samples)
			String mutationKey = new StringBuffer().append(mutationId).append(donorId).append(specimenId)
					.append(sampleId).append(chr).append(start).append(end).append(refBase).append(mutBase).toString();

			Donor donor = donors.get(donorId);

			/*-
			 * Skip this mutation if:
			 * 
			 * 1. Gene affected column is blank (e.g. intergenic mutation)
			 * 2. Gene affected is not in GRch37 (hg19) Ensembl assembly
			 * 3. Gene is NOT protein coding
			 * 4. We threw out this specimen previously (see criteria in readSpecimens)
			 * 5. Not WGS
			 * 6. Mutation is NOT SBS
			 */
			if (StringUtils.isBlank(geneAffected) || !gm.geneInAssembly(geneAffected)
					|| gm.getGene(geneAffected, GeneNaming.ID).isNotProteinCoding()
					|| !donor.containsSpecimen(specimenId) || !Constants.WGS.equals(sequencingStrategy)
					|| !Constants.SBS.equals(mutationType)) {
				continue;
			}

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

				Mutation newMutation = new Mutation(mutationId, donorId, specimenId, sampleId, matchedSampleId, chr,
						start, end, refBase, mutBase, totalReadCount, mutantAlleleReadCount, sequencingStrategy,
						consequenceType, geneAffected);

				newMutation.addRawLine(StringUtils.join(tsvRecord.toList(), '\t'));

				if (Constants.SBS.equals(mutationType)) {
					TriSeq triSeq = ge.getTrinucleotideContext(chr, start);
					newMutation.setTriSeq(triSeq);
				}

				// Retrieve promoter ranges in null-safe manner in case of unusual chr (e.g.
				// "MT")
				List<DnaRange> promoterRanges = MapUtils.getObject(promoterRegionsByChr, chr, new ArrayList<>());
				for (DnaRange promoterRange : promoterRanges) {
					if (promoterRange.isInRange(start, end)) {
						newMutation.setInPromoterRegion(true);

						// No need to keep checking the rest of the promoter ranges
						break;
					}
				}

				// Retrieve CFS ranges in null-safe manner in case of unusual chr (e.g. "MT")
				List<DnaRange> cfsRanges = MapUtils.getObject(cfsRegionsByChr, chr, new ArrayList<>());
				for (DnaRange cfsRange : cfsRanges) {
					if (cfsRange.isInRange(start, end)) {
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

	private void crunchNumbers(String tumorType, String tumorTypeDir, Path outputTimestampFolderPath, Map<String, Donor> donors)
			throws IOException {

		final int highCutoff = 75;
		final int lowCutoff = 25;
		final int zeroCutoff = 0;
		
		Path histogramPath = outputTimestampFolderPath.resolve(tumorType + "_histogram.tsv");
		Path promoterMutationPath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations.tsv");
		Path promoterMutationUniquePath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_unique.tsv");
		
		Path promoterMutationHighPath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_high.tsv");
		Path promoterMutationHighUniquePath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_unique_high.tsv");
		
		Path promoterMutationMidPath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_mid.tsv");
		Path promoterMutationMidUniquePath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_unique_mid.tsv");
		
		Path promoterMutationLowPath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_low.tsv");
		Path promoterMutationLowUniquePath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_unique_low.tsv");
		
		Path promoterMutationZeroPath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_zero.tsv");
		Path promoterMutationZeroUniquePath = outputTimestampFolderPath.resolve(tumorType + "_promoter_mutations_unique_zero.tsv");
		
		Path nonPromoterMutationPath = outputTimestampFolderPath.resolve(tumorType + "_non_promoter_mutations.tsv");
		Path nonPromoterMutationUniquePath = outputTimestampFolderPath.resolve(tumorType + "_non_promoter_mutations_unique.tsv");
		
		Path cfsMutationsPath = outputTimestampFolderPath.resolve(tumorType + "_cfs_mutations.tsv");
		Path cfsMutationsUniquePath = outputTimestampFolderPath.resolve(tumorType + "_cfs_mutations_unique.tsv");
		
		Files.createFile(histogramPath);
		try (BufferedWriter histogramBw = Files.newBufferedWriter(histogramPath);
				BufferedWriter promMutsBw = Files.newBufferedWriter(promoterMutationPath);
				BufferedWriter promMutsUniBw = Files.newBufferedWriter(promoterMutationUniquePath);
				BufferedWriter promMutsHighBw = Files.newBufferedWriter(promoterMutationHighPath);
				BufferedWriter promMutsHighUniBw = Files.newBufferedWriter(promoterMutationHighUniquePath);
				BufferedWriter promMutsMidBw = Files.newBufferedWriter(promoterMutationMidPath);
				BufferedWriter promMutsMidUniBw = Files.newBufferedWriter(promoterMutationMidUniquePath);
				BufferedWriter promMutsLowBw = Files.newBufferedWriter(promoterMutationLowPath);
				BufferedWriter promMutsLowUniBw = Files.newBufferedWriter(promoterMutationLowUniquePath);
				BufferedWriter promMutsZeroBw = Files.newBufferedWriter(promoterMutationZeroPath);
				BufferedWriter promMutsZeroUniBw = Files.newBufferedWriter(promoterMutationZeroUniquePath);
				BufferedWriter nonPromMutsBw = Files.newBufferedWriter(nonPromoterMutationPath);
				BufferedWriter nonPromMutsUniBw = Files.newBufferedWriter(nonPromoterMutationUniquePath);
				BufferedWriter cfsMutsBw = Files.newBufferedWriter(cfsMutationsPath);
				BufferedWriter cfsMutsUniBw = Files.newBufferedWriter(cfsMutationsUniquePath)) {

			histogramBw.write(
					"donorId" + "\t" + "Num specimens" + "\t" + "Num samples" + "\t" + "Num specimens with mutation"
							+ "\t" + "Num specimens with expression" + "\t" + "Num specimens with both" + "\t"
							+ "Num muts" + "\t" + "Num promMuts" + "\t" + "Num promMuts with expression >=" + highCutoff
							+ "th percentile" + "\t" + "Num promMuts with expression >" + lowCutoff + " and <"
							+ highCutoff + "\t" + "Num promMuts with expression <=" + lowCutoff + "\t"
							+ "Num promMuts with expression =" + zeroCutoff + "\t" + "Num CFS mutations");
			histogramBw.newLine();

			// Grab mutation file header
			CSVParser csvParser = openMutationFile(tumorTypeDir);
			List<String> headerElements = new ArrayList<>(csvParser.getHeaderNames());
			csvParser.close();
			headerElements.add("Ref_Tri");
			String header = StringUtils.join(headerElements, '\t');
			promMutsBw.write(header);
			promMutsBw.newLine();
			promMutsHighBw.write(header);
			promMutsHighBw.newLine();
			promMutsHighUniBw.write(header);
			promMutsHighUniBw.newLine();
			promMutsMidBw.write(header);
			promMutsMidBw.newLine();
			promMutsLowBw.write(header);
			promMutsLowBw.newLine();
			promMutsZeroBw.write(header);
			promMutsZeroBw.newLine();
			nonPromMutsBw.write(header);
			nonPromMutsBw.newLine();
			cfsMutsBw.write(header);
			cfsMutsBw.newLine();
			cfsMutsUniBw.write(header);
			cfsMutsUniBw.newLine();

			for (Map.Entry<String, Donor> entry : donors.entrySet()) {
				String donorId = entry.getKey();
				Donor donor = entry.getValue();

				// Only work with donors that have a specimen that has BOTH (based on discussion
				// with Raoof on 5/8/2022)
				if (!donor.hasMutationAndExpressionData()) {
					continue;
				}

				histogramBw.write(donorId);
				histogramBw.write("\t" + donor.getNumSpecimens());
				histogramBw.write("\t" + donor.getNumSamples());
				histogramBw.write("\t" + donor.getNumSpecimensWithMutationData());
				histogramBw.write("\t" + donor.getNumSpecimensWithExpressionData());
				histogramBw.write("\t" + donor.getNumSpecimensWithBoth());

				Map<String, Mutation> mutations = donor.getMutations();
				Set<Mutation> promoterMutations = donor.getPromoterMutations();
				Set<Mutation> nonPromoterMutations = donor.getNonPromoterMutations();

				// Get promoter mutations in high expressed genes (highCutoff <= HERE)
				Set<Mutation> promoterMutationsInHighExpressedGenes = donor.getPromoterMutationsInExpressedGenes(highCutoff, Operator.GREATERTHANOREQUAL);

				// Get promoter mutations in mid expressed genes (lowCutoff < HERE < highCutoff)
				Set<Mutation> promoterMutationsInMidExpressedGenes = donor.getPromoterMutationsInExpressedGenes(highCutoff, Operator.LESSTHAN);
				promoterMutationsInMidExpressedGenes.addAll(donor.getPromoterMutationsInExpressedGenes(lowCutoff, Operator.GREATERTHAN));

				// Get promoter mutations in low expressed genes (HERE <= lowCutoff)
				Set<Mutation> promoterMutationsInLowExpressedGenes = donor.getPromoterMutationsInExpressedGenes(lowCutoff, Operator.LESSTHANOREQUAL);

				// Get promoter mutations in NON-expressed genes (HERE <= zeroCutoff)
				Set<Mutation> promoterMutationsInZeroExpressedGenes = donor.getPromoterMutationsInExpressedGenes(zeroCutoff, Operator.LESSTHANOREQUAL);
				
				Set<Mutation> cfsMutations = donor.getCfsMutations();

				histogramBw.write("\t" + mutations.size());
				histogramBw.write("\t" + promoterMutations.size());
				histogramBw.write("\t" + promoterMutationsInHighExpressedGenes.size());
				histogramBw.write("\t" + promoterMutationsInMidExpressedGenes.size());
				histogramBw.write("\t" + promoterMutationsInLowExpressedGenes.size());
				histogramBw.write("\t" + promoterMutationsInZeroExpressedGenes.size());
				histogramBw.write("\t" + cfsMutations.size());
				histogramBw.newLine();

				writeMutationsToFiles(promoterMutations, promMutsBw, promMutsUniBw);

				writeMutationsToFiles(promoterMutationsInHighExpressedGenes, promMutsHighBw, promMutsHighUniBw);
				
				writeMutationsToFiles(promoterMutationsInMidExpressedGenes, promMutsMidBw, promMutsMidUniBw);
				
				writeMutationsToFiles(promoterMutationsInLowExpressedGenes, promMutsLowBw, promMutsLowUniBw);
				
				writeMutationsToFiles(promoterMutationsInZeroExpressedGenes, promMutsZeroBw, promMutsZeroUniBw);
				
				writeMutationsToFiles(nonPromoterMutations, nonPromMutsBw, nonPromMutsUniBw);
				
				writeMutationsToFiles(cfsMutations, cfsMutsBw, cfsMutsUniBw);

			}
		}
	}

	private void writeMutationsToFiles(Set<Mutation> promoterMutations, BufferedWriter allWriter, BufferedWriter uniqueWriter)
			throws IOException {
		for (Mutation mutation : promoterMutations) {
			String str = mutation.getRawLines().get(0);
			String triSeqWithMut = mutation.getTriSeqWithMut();
			uniqueWriter.write(str + "\t" + triSeqWithMut);
			uniqueWriter.newLine();
			for (String line : mutation.getRawLines()) {
				allWriter.write(line + "\t" + triSeqWithMut);
				allWriter.newLine();
			}
		}
	}

	private CSVParser openMutationFile(String tumorTypeDir) throws IOException {
		Path icgcDir = Paths.get(Constants.WORKING_ICGC_DIR);
		File mutationsFile = icgcDir.resolve(tumorTypeDir).resolve("simple_somatic_mutation.open.tsv.gz").toFile();

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
		// Read in all fragile sites (from Ma et al. 2012 Int J Mol Sci)
		// Cytoband locations from UCSC cytoband.txt.gz hg19

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
