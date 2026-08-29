# Reproducibility Guide

## Project root

All commands below assume the repository has this structure:

```text
TP53-DNA-Sequence-Analyzer/
├── Code/
├── Data/
├── Documentation/
├── Plots/
└── Results/
```

The scripts resolve the repository root from their own location, so they do not require the terminal to start in the repository directory.

## Python

Install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Run:

```bash
python Code/Python/tp53_analysis.py
```

The script regenerates:

```text
Results/Python/TP53_results.txt
Plots/Python/TP53_nucleotide_composition.png
Plots/Python/TP53_gc_at_content.png
```

The Python workflow was tested from outside the repository working directory to verify that it does not depend on `cwd`.

## Linux

Run:

```bash
bash Code/Linux/tp53_analysis.sh
```

The script regenerates deterministic local checks under:

```text
Results/Linux/Reproduced/
├── TP53_nucleotide_composition.txt
├── TP53_physicochemical.txt
├── TP53_reference_validation.txt
└── run_manifest.txt
```

It deliberately does **not** overwrite the historical files directly under `Results/Linux/`. This distinction is important: those files contain the original Linux analysis, including recorded BLASTP results and outputs from external bioinformatics tools.

`run_manifest.txt` records the Python/Biopython versions used for the reproduction run and SHA256 checksums for the four core input FASTA files and the reproduced outputs.

### BLASTP

BLASTP is not treated as a deterministic local dependency because its result depends on the external database contents and query date. The repository therefore preserves the original BLASTP result as recorded evidence. A future fresh BLASTP run should be reported separately with its database/version/date rather than silently replacing the historical result.

## R / RStudio

Install the packages listed in `R_packages.txt`.

Run:

```bash
Rscript Code/R/DNA_analysis.R
```

or open the script in RStudio.

The R workflow resolves `Data/`, `Results/R/`, and `Plots/R/` from the repository root and regenerates the R result file and six static plots.

## Reference validation

The project uses:

- NCBI transcript: `NM_000546.6`
- NCBI protein: `NP_000537.3`
- UniProt: `P04637 (P53_HUMAN)`

The supplied CDS is 1182 nt and corresponds to transcript coordinates 143–1324. The translated CDS and supplied protein reference are 393 aa and match the supplied UniProt reference sequence exactly.

## Historical versus reproduced outputs

The repository intentionally distinguishes between:

1. **Historical recorded outputs** — preserved from the original Linux workflow, including external BLASTP results.
2. **Reproduced deterministic outputs** — generated locally by the cleaned scripts.

This prevents a reproducibility script from destroying evidence of the original analysis while still providing a clean way to verify the core calculations.
