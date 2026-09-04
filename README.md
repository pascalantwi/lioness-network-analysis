# LIONESS Network Analysis

Gene expression network project with Dr Tyler Grines. We are building this pipeline using osteosarcoma data as practice before moving to lupus. We are rebuilding patient specific networks from raw GEO data using LIONESS with Pearson correlation.

## Background

LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples) estimates an individual network for each patient in a dataset, rather than a single average network for the whole population. We are using the osteosarcoma dataset from Kuijjer et al. (2019, BMC Cancer) as a proof of concept, since it is well documented and the original authors' pre packaged results are available for comparison. Once this pipeline is validated, we will apply the same approach to lupus gene expression data.

## Data

- Source: GEO accessions GSE42352 (superseries) and GSE33382 (subseries containing the 84 high grade osteosarcoma pre chemotherapy biopsies)
- Platform: GPL10295, Illumina human-6 v2.0 expression beadchip
- Outcome variable: metastasis within 5 years (yes/no), 53 patients with known follow up
- We reconstruct gene expression data from the raw GEO series matrix rather than using the pre packaged `OSdata.RData` file directly, so the pipeline can be adapted for datasets that are not already processed, such as the lupus data.

## What this repository contains

- `RA.R` - full pipeline, from raw GEO data through LIONESS, LIMMA differential co-expression, and network visualization. Includes both our own raw data reconstruction and a parallel run on the official `OSdata.RData` for validation.

## Pipeline overview

1. Load raw GEO series matrix and platform annotation
2. Filter to the relevant patient samples and map probes to gene symbols
3. Collapse duplicate gene symbols, keeping the highest variance probe
4. Subset to the most variable genes
5. Run LIONESS to estimate single sample networks
6. Use LIMMA to find edges that differ significantly between patient outcome groups
7. Visualize the top differential edges as a network graph

## Requirements

R packages: `GEOquery`, `AnnotationDbi`, `lionessR`, `igraph`, `reshape2`, `limma`, `pheatmap`

## Status

Building. Osteosarcoma proof of concept is running and being validated against the original authors' data.
