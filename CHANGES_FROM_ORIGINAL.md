# Changes from MattStata/RNAseqWebinar

This rewrite preserves the original workflow stages while removing scheduler-specific execution.

| Original file | Replacement | Main changes |
|---|---|---|
| `01_run_fastqc_multiqc.slurm` | `01_run_fastqc_multiqc.sh` | Removes `#SBATCH`/modules; uses Conda PATH and GNU Parallel with configurable local concurrency. |
| `02_run_trimmomatic.slurm` | `02_run_trimmomatic.sh` | Removes SLURM; reads a manifest; runs Trimmomatic PE or SE per sample; emits a trimmed manifest. |
| `03_run_hisat2.slurm` | `03_run_hisat2.sh` | Removes SLURM; supports HISAT2 `-1/-2` for PE and `-U` for SE; local index reuse and BAM sorting/indexing. |
| `04_run_salmon.slurm` | `04_run_salmon.sh` | Removes SLURM; supports Salmon `-1/-2` for PE and `-r` for SE; supports supplied transcript references or genome+GFF3 reference preparation. |
| `QuantifyAllBAM.py` | `QuantifyAllBAM.py` | Adds manifest-aware mixed SE/PE BAM counting while retaining `--single` compatibility. |
| `QuantifySalmon.py` | `QuantifySalmon.py` | Self-contained matrix aggregation from Salmon `quant.sf` / `quant.genes.sf`. |
| `VisualizeResults.py` | `VisualizeResults.py` | Self-contained PCA/correlation/t-SNE plotting with the same general downstream role. |

Additional files:

- `00_make_manifest.py`: infers SE/PE layout from common FASTQ names.
- `config.sh`: paths, references, and workstation resources.
- `config_dgx_spark.sh`: DGX Spark 20-core concurrency preset.
- `environment.yml`: regular Conda environment.
- `environment_dgx_spark.yml`: native `linux-aarch64` DGX Spark environment.
- `prepare_salmon_reference.py`: uses gffread for strand-aware transcript extraction and creates transcript→gene mapping from genome+GFF3/GTF when needed.
- `validate_environment.sh`: command/import validation.
- `run_all.sh`: optional sequential driver.
