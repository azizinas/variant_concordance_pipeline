# Concordia - A variant concordance pipeline
A variant concordance pipeline for benchmarking NGS callers against truth sets.
Concordia automates the comparison of VCF files using RTG Tools' vcfeval, annotates discordant variants with likely failure modes, and produces an interactive HTML report.

## Background
Evaluating variant caller performance requires more than counting matches and mismatches, discordant variants need context. Are they in low-complexity regions? High GC content? Caller-filtered? Concordia answers these questions systematically, turning a raw vcfeval result into an annotated, interpretable report. The pipeline supports multiple samples via a CSV samplesheet and processes them in parallel, making it practical for batch benchmarking runs.

## Pipeline overview
Concordia is structured as three Python modules orchestrated by a Nextflow workflow:

main.nf — Nextflow workflow entry point, parameter handling
Reference prep & vcfeval — FASTA indexing, RTG vcfeval execution
Discordant annotation — failure mode annotation per variant
Report prep — TSV generation, Quarto HTML report

### Failure mode annotations

Each discordant variant is annotated with one or more of the following flags:

| Flag | Description |
|------|-------------|
| `GT_DISCORDANCE` | Variant called but genotype differs from truth |
| `HIGH_GC` | Variant falls in a high GC/AT content region (via `hg38.gc5Base.bw`) |
| `LOW_COMPLEXITY` | Variant overlaps a low-complexity or repeat region |
| `SEGDUP` | Variant falls within a segmental duplication |
| `CALLER_FILTER` | Variant was soft-filtered by the caller |
| `REGION_BASED` | Variant falls outside callable regions |

Annotations are computed using `pyBigWig` for GC content and `pybedtools` for interval-based overlap.

## Inputs

samplesheet (CSV, required): columns — sample_name, run_id, baseline, query
targets_bed (BED, optional): target regions file — required for WES, targeted panel, and amplicon assays
mode: germline or somatic
build: hg19 or hg38
assay_type: wgs, wes, targeted_panel, or amplicon

## Outputs

vcfeval/ — raw vcfeval output including summary.txt
discordance_modes/ — TSV files per failure mode
report_inputs/ — metadata, counts, and performance TSVs
Report.html — interactive concordance report
Reference Data

Genome FASTA downloaded at runtime from UCSC (hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz, hg38: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz)
GC content: gc5Base.bw from UCSC Genome Browser track (hg19/hg38)
Low complexity regions and segmental duplications: GIAB Genome Stratifications v3.6 (https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/)

## Container

Base image: realtimegenomics/rtg-tools:3.13
Python dependencies: requests, pandas, cyvcf2, pybedtools, pyBigWig
R dependencies: DT, UpSetR
Additional tools: bcftools, bedtools, Quarto
