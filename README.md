# *Salmonella* Enteritidis genomic epidemiology (Saudi Arabia + global context)

Code to reproduce the main figures of a *Salmonella* Enteritidis population-genomics
study combining a Saudi Arabian collection with related publicly available genomes.

## Files

**`scripts/`**

| File | Description |
|------|-------------|
| `Enteritidis_plot.R` | Full R pipeline producing Figures 1–5 |

**`Figure1/`** (input data read by the script)

| File | Description |
|------|-------------|
| `Enteritidis_snp_sites.aln` | Core-genome SNP alignment (397 genomes) |
| `metadata_baps.tsv` | Per-genome metadata: hierBAPS group, source, country, city, date |
| `amrfinder_all.tsv` | AMRFinderPlus results (all genomes) |
| `plasmid_freq.csv` | Per-genome plasmid reference coverage (%) |
| `plasmid_replicon_by_contig.csv` | Replicon type of each representative plasmid |
| `prophage_freq.tsv` | ST64B prophage coverage per genome |
| `BAPS_Age_clock_results.xlsx` | BEAST MRCA age and clock-rate estimates |
| `baps{1,2,3,6,7}_country_beast1.tree` | Time-scaled BEAST trees |
| `baps{1,2,3,4,6,7,8}.filtered_polymorphic_sites.fasta` | Per-group SNP alignments |

## Running

Run from the repository root so the `Figure1/` input paths resolve:

```bash
Rscript scripts/Enteritidis_plot.R
```

Figures are written as SVG into `Figure1/`.

---
## Contact

For questions about this repository, contact Ge Zhou (ge.zhou@kaust.edu.sa; PhD
candidate) at the [Infectious Disease Epidemiology Lab](https://ide.kaust.edu.sa/), KAUST.



