# Health-associated gut bacteriocins target TLR4 to suppress intestinal inflammation

Scripts for mapping class II bacteriocin sequences to metagenomic and metatranscriptomic reads, and for differential abundance analysis across cohorts.

## Contents

```
Multi-omics_mapping/
├── mapping.sh                          # PALADIN mapping pipeline
├── get_count_table.py                  # Build a count table from samtools idxstats
├── analysis_code.R                  # Diversity + DESeq2 differential abundance analysis
├── 00.data/
│   ├── metadata.csv                    # Sample metadata (dataset, sample, group)
│   └── overall_High.id_0.95_co_0.95.fas  # Clustered class II bacteriocin sequences
└── 00.output/
     └── counts.txt                      # Combined per-sample count table
```

## Requirements

- PALADIN
- samtools
- Python 3
- R (>= 4.0) with the packages `dplyr`, `DESeq2`, and `tibble`

## Input files

- `00.data/overall_High.id_0.95_co_0.95.fas`: clustered class II bacteriocin nucleotide sequences used as the PALADIN reference.
- Clean paired-end reads, organized in per-sample directories named `SRR*`.
- `00.data/metadata.csv`: sample metadata with the columns `Datasets`, `Samples`, and `Group` (`Healthy`, `UC`, or `CRC`).

## Pipeline
mapping.sh 

### Step 1: Build the PALADIN index

### Step 2: Map clean reads

### Step 3: Build the count table

### Step 4: Differential abundance analysis
