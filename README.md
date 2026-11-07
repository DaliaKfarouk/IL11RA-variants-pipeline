# IL11RA variants: an independent re-annotation of published disease variants

**Status (4 November 2026):** Modules 0, 1 and 2 are done and run end to end. Modules 3 and 4 are planned and not started. Open data issues are listed in [`docs/deviations_log.md`](docs/deviations_log.md).

## What this project does

People who are born without a working IL-11 receptor (the gene is *IL11RA*) can have fused skull sutures (craniosynostosis) and tooth problems. A systematic review catalogues every published *IL11RA* variant. This repository takes that catalogue and asks one question: **can an independent, scripted workflow reproduce what the papers reported, and what would a standard rare-disease filter miss?**

It generates no new patient data. All inputs are published or public. The aim is a transparent check, with every choice written down.

## Questions answered so far

| # | Question | Module |
|---|---|---|
| Q1 | Does an independent annotation workflow reproduce the published protein change of each variant? Which published disease variants would a standard "rare, protein-altering, high-score" filter pass or miss? How do the papers, ClinVar, gnomAD and AlphaMissense agree? | 1 |
| Q2 | Where do the variants fall in the receptor (signal peptide, domains, WSXWS motif)? | 2 |

Modules 3 (public RNA-seq of IL11 induction and tissue expression) and 4 (druggability and target-safety evidence) are intended to be finished by the end of December and added here 

## Main results so far

All numbers come from `results/` and were produced by the scripts in `R/`.

- **Input:** 31 *IL11RA* variant rows, from 26 families in 12 source reports. One further row (a different gene) is kept apart in `results/m1_other_genes.csv`. 30 rows validated in VariantValidator; 1 did not (VAR11).
- **Reproducing the published protein change (Q1a):** 25 variants could be predicted by `predictCoding`. 21 matched the published name (84.0%). 4 did not (VAR10, VAR11, VAR27, VAR31). These are not yet resolved. I describe each one in `notes/m1_interpretation.md`.
- **Reference residue check:** for 28 of the 29 rows that have a protein name, the amino acid at the published position matches UniProt Q14626. The one exception is VAR11, which is off by one residue.
- **Detectability (Q1b):** the permissive filter F1 and the inclusive filter F3 pass all 31 variants. The standard filter F2 passes 19 of 31 (61.3%, exact 95% CI 42.2% to 78.2%). All 12 variants that F2 misses fail on variant class or AlphaMissense score. None fail on frequency. The misses are 8 missense variants scored "ambiguous" or "likely benign", 3 in-frame deletions or duplications, and 1 splice-region variant.
- **Source agreement:** 12 of 31 variants are in ClinVar. 19 missense variants are scored by AlphaMissense (11 likely pathogenic, 4 ambiguous, 4 likely benign).
- **Protein map (Q2):** 29 variants have a protein position. 18 of them fall in the third fibronectin-like domain (D3, 98 of 422 residues). This is a description of a small set, not a test.

## Repository layout

```
R/                     numbered scripts, run in order (see run_all.R)
data/                  frozen SR table, VariantValidator coordinates, domain table
data/raw/              downloads from ClinVar, gnomAD, AlphaMissense (not committed)
results/               tables written by the scripts
figures/               figures written by the scripts
docs/                  analysis plan, data dictionary, deviations log, project walkthrough
notes/                 dated working log and the written interpretation of each module
```

## How to reproduce

1. Install R (this project was run with R 4.6.1). Open the folder as an RStudio project so relative paths work.
2. Put these files under `data/`:
   - `data/variants_input.csv`: the frozen systematic-review variant table. Never edit it by hand.
   - `data/variants_coords.csv`: GRCh38 coordinates and MANE names from VariantValidator.
   - `data/raw/`: the ClinVar table download, the gnomAD variant export for *IL11RA* (gene ENSG00000137070), and the AlphaMissense search export for UniProt Q14626 (leave it zipped).
   - `data/domains_IL11RA.csv`: start and end of each domain, copied from UniProt Q14626.
3. Run `source("R/00_setup.R")`. It installs any missing packages and builds the folders. The first run is slow because the hg38 genome package is large.
4. Run `source("R/run_all.R")`, or run scripts 01 to 05 one at a time.

Package versions are locked in `renv.lock`. 

## Data sources

| Resource | Version or date | Note |
|---|---|---|
| Systematic-review variant table | frozen 2026-11-04 | checksum in `results/m1_input_freeze.csv` |
| VariantValidator | | GRCh38 coordinates, MANE names |
| ClinVar | downloaded 2026-10-30 | table download |
| gnomAD | v4.1.2 as written in the export file name; exported 2026-10-31 | genome build GRCh38 |
| AlphaMissense | search export for Q14626, file dated 2026-10-30 | licence to be checked before sharing the raw file |
| UniProt Q14626 | | 422 amino acids |
| Transcript models | `TxDb.Hsapiens.UCSC.hg38.knownGene` | see `renv.lock` for the version |

Every file has its checksum in `results/m1_external_sources.csv`. Raw downloads are not committed. Download them again from the sources above.

## Limits you should know about

- The human sample is small and biased towards striking cases. It comes from published reports only.
- I **reproduced** published information with an independent workflow. That is not independent **validation**, and I use the word "reproduced" on purpose.
- AlphaMissense, ClinVar and gnomAD are supporting evidence. They do not prove that a variant causes disease.
- Variants absent from the gnomAD export count as rare in the filters. For three variants (VAR05, VAR17, VAR20) no genomic position exists, so "absent" really means "could not be looked up".
- Four variants did not reproduce, and three SR rows need checking against the papers (VAR11, VAR21, VAR31). Details are in `docs/deviations_log.md`. Numbers in this README will change if the SR rows are corrected and the pipeline is re-run.
- The functional-evidence classes in the protein map come from keyword rules applied to free text. "Other" does not mean "no evidence".
- One analyst ran everything. A second person has not yet re-run it from scratch.

## Relationship to the systematic review

The systematic review is registered on OSF. This repository takes the frozen variant table as input and sends corrections back to the review. It never edits the review's columns. Variant classes follow the review protocol.

## AI assistance

I used Claude to assist with the errors R scripts returned. I ran the scripts, fixed the errors that came up, and checked the outputs against the source files. I am responsible for the content. 

## contact

Dalia Khaled Farouk Kamel 

