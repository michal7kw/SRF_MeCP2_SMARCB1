# Pipeline Execution

1. `run_{sample}.sh`
   - If failed: `continue_{sample}.sh`
2. `run_smarcb1_analysis`
   - `run_smarcb1_analysis_{control}_{sample}_{cell_type}_{metric}.sh`
      - `submit_selected.sh`


`all_mecp2_targets_1.csv` --> joined: 
- `enriched_down_regulated.csv`
- `enriched_not_disregulated.csv`
- `enriched_up_regulated.csv`
(extracted from .xlsx files from email)

`all_mecp2_targets_2.csv` --> subseted (for: both and exon_only): 
- `complete_peak_annotation.csv`
(originated from  `SRF_MeCP2_CUTandTAG/iterative_alternative/results/no_dedup/cpg_enrichment/NSC/broad/cpg_enrichment_2_rep_in_peaks/`)