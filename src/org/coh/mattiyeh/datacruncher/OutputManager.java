package org.coh.mattiyeh.datacruncher;

import java.util.ArrayList;
import java.util.List;

public class OutputManager {
	
	private List<String> metadataHeaderItems;
	
	public OutputManager(int highCutoff, int lowCutoff) {
		metadataHeaderItems = new ArrayList<>();
		metadataHeaderItems.add("donor_id");
		metadataHeaderItems.add("specimens");
		metadataHeaderItems.add("samples");
		metadataHeaderItems.add("spec_w_muts");
		metadataHeaderItems.add("spec_w_exp");
		metadataHeaderItems.add("spec_w_both");
		
		metadataHeaderItems.add("all_muts");
		
		metadataHeaderItems.add("non_prom_muts");
		metadataHeaderItems.add("prom_muts");
		
		metadataHeaderItems.add("non_prom_sbs_muts");
		metadataHeaderItems.add("prom_sbs_muts");
		metadataHeaderItems.add("prom_sbs_muts_w_exp_gte_" + highCutoff);
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
	}
	
	public List<String> getMetadataHeaderItems() {
		return metadataHeaderItems;
	}

}
