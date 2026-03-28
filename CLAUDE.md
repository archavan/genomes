# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Research project for building processed genomic reference databases across multiple species, primarily for comparative genomics. Species include *Styela clava* (scla), *Mus musculus* (mmus), *Ciona robusta* (crob), *Brachionus floridanus* (bflo), and others. Mouse serves as the reference species for functional annotations.

## Running Code

This is an RStudio project (`genomes.Rproj`) — no build system. Scripts are meant to be sourced in R or run from Positron (primary IDE). Quarto documents (`.qmd`) are rendered via Positron or `quarto render`.

There are no tests, no linter configuration, and no CI/CD. Quality assurance is done through Quarto reports with inline validation.

## Code Architecture

### Directory Layout

```
species/<abbrev>/<source-version>/
    code/          # numbered processing scripts
    data/          # raw downloads (gitignored)
    results/       # processed outputs (gitignored)

orthofinder/<date>/
    code/          # comparative analysis scripts
    data/          # species.tsv and genome URLs
    results/       # OrthoFinder outputs + processed tables (gitignored)
```

### Species-Level Processing Pipeline

Scripts are numbered to indicate execution order (e.g., `00a_`, `01a_`, `02a_`). Each species folder follows the same pattern:

1. **Download** — shell scripts fetch genome files from RefSeq/Ensembl/Ghost
2. **GTF parsing** — `rtracklayer` parses annotation files; gene IDs are fixed (LOC/GeneID prefixes resolved, pseudogenes removed)
3. **Gene info extraction** — coordinates, symbols, cross-database IDs consolidated into TSVs
4. **Functional annotation** — eggNOG MAPPER output parsed for GO terms and KEGG IDs; KEGG REST API queried for pathways/modules
5. **Longest proteins** — one protein per gene extracted for OrthoFinder input
6. **Display names** — `11a_generate-gene-display-names.R` creates human-readable labels, incorporating OrthoFinder HOG labels when available
7. **Gene sets** — transcription factors (`12a_`), functional modules (`13a_`)

### OrthoFinder Comparative Analysis (`orthofinder/<date>/`)

```
00b_prep-input-proteomes.R      → standardize FASTA headers
01a_run-orthofinder.sh          → execute OrthoFinder binary
02a_process-output-files.R      → parse orthogroups, extract pairwise orthologs
03a_plot-gene-trees.R           → visualize gene trees
04a_sanity-checks.qmd           → validation report
05a_extract-pairwise-orthologs.R
06a_write-hog-labels_<sp>.R     → per-species HOG label generation (utils in utils_hog-labels.R)
```

HOG labels from OrthoFinder are fed back into each species' display name generation.

### Key Output Tables

All primary outputs are TSV files. Standard columns:
- `gene_id` — stable identifier (prefixed: `ncbi.`, `kegg.`, etc.)
- `gene_symbol` — human-readable name
- `orthogroup`, `hog_label` — orthology assignments from OrthoFinder
- Cross-database IDs: `kegg_gene_id`, `ko_id`, `uniprot_id`
- Annotations: `go_terms`, `kegg_pathways`, `description`

### Formatting

R code is formatted automatically by [Air](https://posit-dev.github.io/air/) on save. Do not manually reformat R code.

### Key R Libraries

`tidyverse`, `rtracklayer`, `Biostrings`, `GenomicRanges`, `GO.db`, `KEGGREST`, `ape`, `phytools`, `here`, `janitor`

### Important Conventions

- `here::here()` is used for all file paths relative to the project root
- Data and results directories are gitignored; only code and small metadata files (e.g., `species.tsv`, curated TSVs like `functional-modules_gut.tsv`) are tracked
- Mouse functional modules (`mmus/gene-modules/`) serve as the reference for cross-species mapping
