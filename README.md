# CORRAL

A Nextflow pipeline that aligns sequencing reads against a marker-gene reference database and summarizes taxon-level abundance from the alignments.

## Overview

This pipeline aligns short reads to a marker-gene reference database with `bowtie2`, then summarizes the resulting alignments into taxon abundance tables using `marker_alignments`. It is used within VEuPathDB's metagenomic/eukaryotic-community profiling workflows to estimate the taxonomic composition of a sample from marker-gene read alignments. Reads can be pulled from SRA by run accession or read directly from local FASTQ files; for each sample, reads are aligned to a reference database, per-alignment statistics are computed, alignments are summarized into per-sample taxon tables, and the per-sample tables are finally combined into one summary matrix.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- Docker (enabled by default via `nextflow.config`)

The pipeline runs on the `veupathdb/corral:latest` container image, which bundles `bowtie2`, `samtools`, the SRA toolkit (`fasterq-dump`), and the `marker_alignments` package.

## Usage

```
nextflow run VEuPathDB/CORRAL \
  -r main \
  -entry sra \
  --inputPath /path/to/runAccessions.tsv \
  --refdb /path/to/database/refdb \
  --markerToTaxonPath /path/to/marker-to-taxon-id.tsv \
  --libraryLayout paired \
  --resultDir /path/to/output \
  -resume
```

The pipeline has two named workflow entry points:

- **`sra`** — reads a TSV of run accessions from `params.inputPath` (a `run_accession` column), downloads each run's reads from SRA with `fasterq-dump`, then aligns, scores, and summarizes them.
- **`local`** — reads paired FASTQ files directly from the directory given in `params.inputPath` (matched as `*_{1,2}.fastq`), then aligns, scores, and summarizes them without a download step.

Both entry points converge on the same alignment/summarization logic: `bowtie2` alignment, per-sample alignment statistics, per-sample taxon summarization, and a final combined summary table.

## Key Parameters

| Parameter | Description |
| --- | --- |
| `params.downloadMethod` | `sra` or `local`; selects whether reads are downloaded from SRA or read from local files (used only when running without an explicit `-entry`). |
| `params.inputPath` | For `sra`: path to a TSV listing `run_accession` values. For `local`: path to a directory of paired FASTQ files. |
| `params.libraryLayout` | `single` or `paired`; selects which `bowtie2` alignment template is used. |
| `params.refdb` | Path prefix passed to `bowtie2 -x` for the marker-gene reference database. |
| `params.bowtie2Command` | The `bowtie2` command line (including flags) used for alignment. |
| `params.alignmentStatsCommand` | Command used to compute alignment statistics (e.g. `samtools stats`). |
| `params.summarizeAlignmentsCommand` | The `marker_alignments` command line (including filter thresholds) used to summarize alignments into taxon calls. |
| `params.markerToTaxonPath` | Path to the marker-to-taxon-ID mapping TSV used by `marker_alignments`. |
| `params.summaryColumn` | Value column used when building the final summary matrix (e.g. `cpm`). |
| `params.summaryFormat` | Output format for the final summary table (e.g. `matrix`). |
| `params.resultDir` | Directory where alignment stats, per-sample taxon summaries, and the final combined table are published. |

## Output

Under `params.resultDir`: an `alignmentStats/` directory with per-sample alignment statistics, a `summarizedAlignments/` directory with per-sample `<sample>.taxa.tsv` taxon summaries, and a combined `<summaryColumn>.<summaryFormat>.tsv` file with taxon abundance across all samples.
