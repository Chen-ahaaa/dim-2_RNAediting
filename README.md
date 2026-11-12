# dim-2_RNAediting

Analysis scripts accompanying research on DIM-2 and RNA editing, covering RNA editing detection, DNA methylation analysis, mutation filtering, and characterization of DIM-2 homologous coding sequences.

## Scripts

| Directory | Script | Description |
| --- | --- | --- |
| `RNA_editing/` | `RNA_editing_pipeline.sh` | Processes RNA-seq reads and detects candidate A-to-I editing sites using REDItools, with optional DNA BAM input. |
| | `anno_reditool_res.py` | Filters editing candidates by gene orientation and substitution type, and annotates sites using SnpEff. |
| | `merge_reditools_re.py` | Merges results across DNA comparisons and strand-specific analyses for each RNA sample. |
| `DNA_methylation/` | `WGBS_analysis.sh` | Performs Bismark alignment, deduplication, methylation extraction, lambda non-conversion estimation, and methylation summaries by genomic region class. |
| | `config.example.sh` | Provides paths and parameters for the WGBS workflow. |
| | `methy_by_window.sh` | Calculates per-sample methylation levels from Bismark CX reports within supplied BED windows. |
| `Mutation_filter/` | `detect_snp.sh` | Processes existing VCF and BAM files, quantifies allele support, and coordinates mutation filtering and downstream inspection. |
| | `filter_spore_mut_v8.py` | Filters and classifies candidate mutations using parental and progeny read support and segregation patterns; supports random spores and tetrads. |
| `Species_tree_with_stop_codon/` | `filter_dim2.py` | Selects representative DIM-2 CDS sequences per species and exports sequences, species names, and terminal stop codon types. |

## Requirements and usage

Dependencies vary by workflow and include Trimmomatic, HISAT2, SAMtools, Picard, REDItools, SnpEff, Bismark, Bowtie2, BEDTools, BCFtools, VCFtools, and VarScan. Python packages include `pandas`, `pysam`, and `pyranges`. The legacy REDItools step requires a separate Python 2.7 environment.

Run shell workflows in Linux/Bash after preparing input data and reference files and updating paths and sample-specific parameters. Some auxiliary scripts referenced by `detect_snp.sh` are not included and must be supplied separately. `filter_dim2.py` reads `fasta.txt` from the working directory and expects `organism_name` metadata in FASTA headers.

For the WGBS workflow, copy and edit the configuration template before running:

```bash
cp DNA_methylation/config.example.sh DNA_methylation/config.sh
# Edit config.sh to specify paths and parameters.
bash DNA_methylation/WGBS_analysis.sh DNA_methylation/config.sh
```
