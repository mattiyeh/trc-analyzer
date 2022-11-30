package org.coh.mattiyeh.datacruncher;

import java.util.ArrayList;
import java.util.List;

public class OutputManager {

	private final List<String> metadataHeaderItems;
	private final List<String> summaryHeaderItems;

	public OutputManager(int highCutoff, int lowCutoff) {
		metadataHeaderItems = new ArrayList<>();
		metadataHeaderItems.add("donor_id");
		metadataHeaderItems.add("specimens");
		metadataHeaderItems.add("samples");
		metadataHeaderItems.add("spec_w_muts");
		metadataHeaderItems.add("spec_w_exp");
		metadataHeaderItems.add("spec_w_both");

		metadataHeaderItems.add("muts");
		metadataHeaderItems.add("sbs_muts");
		metadataHeaderItems.add("sbs_muts_C_to_T");
		metadataHeaderItems.add("sbs_muts_C_to_A");
		metadataHeaderItems.add("sbs_muts_C_to_G");
		metadataHeaderItems.add("sbs_muts_T_to_C");
		metadataHeaderItems.add("sbs_muts_T_to_A");
		metadataHeaderItems.add("sbs_muts_T_to_G");
		metadataHeaderItems.add("indel_muts");
		metadataHeaderItems.add("mbs_muts");

		metadataHeaderItems.add("non_prom_muts");
		metadataHeaderItems.add("prom_muts");

		metadataHeaderItems.add("non_prom_sbs_muts");
		metadataHeaderItems.add("prom_sbs_muts");
		metadataHeaderItems.add("prom_sbs_muts_C_to_T");
		metadataHeaderItems.add("prom_sbs_muts_C_to_A");
		metadataHeaderItems.add("prom_sbs_muts_C_to_G");
		metadataHeaderItems.add("prom_sbs_muts_T_to_C");
		metadataHeaderItems.add("prom_sbs_muts_T_to_A");
		metadataHeaderItems.add("prom_sbs_muts_T_to_G");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff);
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff + "_C_to_T");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff + "_C_to_A");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff + "_C_to_G");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff + "_T_to_C");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff + "_T_to_A");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff + "_T_to_G");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
		metadataHeaderItems.add("prom_sbs_muts_w_exp_lte_" + lowCutoff);
		metadataHeaderItems.add("prom_sbs_muts_w_zero");

		metadataHeaderItems.add("non_prom_indel_muts");
		metadataHeaderItems.add("prom_indel_muts");
		metadataHeaderItems.add("prom_indel_muts_w_exp_gte_" + highCutoff);
		metadataHeaderItems.add("prom_indel_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
		metadataHeaderItems.add("prom_indel_muts_w_exp_lte_" + lowCutoff);
		metadataHeaderItems.add("prom_indel_muts_w_zero");

		metadataHeaderItems.add("non_prom_mbs_muts");
		metadataHeaderItems.add("prom_mbs_muts");
		metadataHeaderItems.add("prom_mbs_muts_w_exp_gte_" + highCutoff);
		metadataHeaderItems.add("prom_mbs_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
		metadataHeaderItems.add("prom_mbs_muts_w_exp_lte_" + lowCutoff);
		metadataHeaderItems.add("prom_mbs_muts_w_zero");

		metadataHeaderItems.add("non_cfs_muts");
		metadataHeaderItems.add("cfs_muts");

		metadataHeaderItems.add("non_cfs_sbs_muts");
		metadataHeaderItems.add("cfs_sbs_muts");

		metadataHeaderItems.add("non_cfs_indel_muts");
		metadataHeaderItems.add("cfs_indel_muts");

		metadataHeaderItems.add("non_cfs_mbs_muts");
		metadataHeaderItems.add("cfs_mbs_muts");

		summaryHeaderItems = new ArrayList<>();
		summaryHeaderItems.add("tumor_type");
		summaryHeaderItems.add("specimens");
		summaryHeaderItems.add("samples");
		summaryHeaderItems.add("spec_w_muts");
		summaryHeaderItems.add("spec_w_exp");
		summaryHeaderItems.add("spec_w_both");

		summaryHeaderItems.add("muts");
		summaryHeaderItems.add("sbs_muts");
		summaryHeaderItems.add("indel_muts");
		summaryHeaderItems.add("mbs_muts");

		summaryHeaderItems.add("non_prom_muts");
		summaryHeaderItems.add("prom_muts");

		summaryHeaderItems.add("non_prom_sbs_muts");
		summaryHeaderItems.add("prom_sbs_muts");
		summaryHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff);
		summaryHeaderItems.add("prom_sbs_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
		summaryHeaderItems.add("prom_sbs_muts_w_exp_lte_" + lowCutoff);
		summaryHeaderItems.add("prom_sbs_muts_w_zero");

		summaryHeaderItems.add("non_prom_indel_muts");
		summaryHeaderItems.add("prom_indel_muts");
		summaryHeaderItems.add("prom_indel_muts_w_exp_gte_" + highCutoff);
		summaryHeaderItems.add("prom_indel_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
		summaryHeaderItems.add("prom_indel_muts_w_exp_lte_" + lowCutoff);
		summaryHeaderItems.add("prom_indel_muts_w_zero");

		summaryHeaderItems.add("non_prom_mbs_muts");
		summaryHeaderItems.add("prom_mbs_muts");
		summaryHeaderItems.add("prom_mbs_muts_w_exp_gte_" + highCutoff);
		summaryHeaderItems.add("prom_mbs_muts_w_exp_gt_" + lowCutoff + "_lt_" + highCutoff);
		summaryHeaderItems.add("prom_mbs_muts_w_exp_lte_" + lowCutoff);
		summaryHeaderItems.add("prom_mbs_muts_w_zero");

		summaryHeaderItems.add("non_cfs_muts");
		summaryHeaderItems.add("cfs_muts");

		summaryHeaderItems.add("non_cfs_sbs_muts");
		summaryHeaderItems.add("cfs_sbs_muts");

		summaryHeaderItems.add("non_cfs_indel_muts");
		summaryHeaderItems.add("cfs_indel_muts");

		summaryHeaderItems.add("non_cfs_mbs_muts");
		summaryHeaderItems.add("cfs_mbs_muts");
	}

	public List<String> getMetadataHeaderItems() {
		return new ArrayList<>(metadataHeaderItems);
	}
	
	public List<String> getSummaryHeaderItems() {
		return new ArrayList<>(summaryHeaderItems);
	}

}
