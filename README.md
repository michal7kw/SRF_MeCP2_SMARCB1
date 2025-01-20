# Pipeline Execution

1. `run_{sample}.sh`
   - If failed: `continue_{sample}.sh`
2. `run_smarcb1_analysis`
   - `run_smarcb1_analysis_{control}_{sample}_{cell_type}_{metric}.sh`
      - `submit_selected.sh`
