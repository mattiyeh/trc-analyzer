package org.coh.mattiyeh.datacruncher;

public class Constants {

	static final String WORKING_DIR = "C:\\Users\\" + ("CONTRA".equals(Launcher.getHostname()) ? "marca\\Google " : "attiyehm\\My ")
			+ "Drive\\City of Hope\\Projects\\TRC signatures\\";
	
	static final String DATA_FILES_DIR = WORKING_DIR + "data_files\\";
	static final String OUTPUT_DIR = WORKING_DIR + "icgc-datasets-outputs\\";
	
	static final String DONORS_FILENAME = "donor.tsv.gz";
	static final String SPECIMENS_FILENAME = "specimen.tsv.gz";
	static final String SAMPLES_FILENAME = "sample.tsv.gz";
	static final String MUTATIONS_FILENAME = "simple_somatic_mutation.open.tsv.gz";

	static final String[] TUMOR_TYPES2 = { "bladder", "bone", "breast", "cervix", "colorectal", "esophagus",
			"gallbladder", "kidney", "liver", "lung", "melanoma", "ovary", "pancreas", "pnet", "prostate", "stomach",
			"uterus" };
	
	static final String[] TUMOR_TYPES = { "cervix" };

	static final String WORKING_ICGC_DIR = WORKING_DIR + "icgc-datasets\\";
	static final String ICGC_SAMPLE_ID = "icgc_sample_id";
	static final String ICGC_DONOR_ID = "icgc_donor_id";
	static final String ICGC_SPECIMEN_ID = "icgc_specimen_id";
	static final String SPECIMEN_TYPE = "specimen_type";
	static final String SPECIMEN_DONOR_TREATMENT_TYPE = "specimen_donor_treatment_type";

	static final String SEQUENCING_STRATEGY = "sequencing_strategy";

	static final String CHEMOTHERAPY = "chemotherapy";
	static final String CHEMORADIATION = "combined chemo+radiation therapy";

	static final String PRIMARY = "Primary tumour";
	//static final String CELL_LINE = "Cell line";
	//static final String NORMAL = "Normal";
	//static final String XENOGRAFT = "Xenograft";
	//static final String RECURRENT = "Recurrent";
	//static final String METASTATIC = "Metastatic";
	
	//static final String BLOOD = "blood";
	static final String LYMPH_NODE = "lymph node";
	static final String SOLID_TISSUE = "solid tissue";
	static final String ADDITIONAL_NEW_PRIMARY = "additional new primary";
	static final String OTHER = "other";

	public static final String SBS = "single base substitution";
	public static final String INSERTION = "insertion of <=200bp";
	public static final String DELETION = "deletion of <=200bp";
	public static final String MBS = "multiple base substitution (>=2bp and <=200bp)";
	
	static final String WGS = "WGS";

	private Constants() {
		throw new IllegalStateException("Utility class");
	}

}
