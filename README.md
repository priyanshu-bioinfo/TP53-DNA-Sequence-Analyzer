# Human TP53 Nucleotide Sequence Validation and Feature Analysis

A reproducible bioinformatics project analyzing the **Homo sapiens TP53** sequence using **Python, Linux, and R/RStudio**.

## Project objective

The project analyzes the human TP53 transcript, its annotated coding sequence (CDS), and the corresponding protein. It combines programmable sequence analysis, command-line bioinformatics, database validation, protein-property analysis, ORF/frame analysis, and visualization.

## Reference records

- NCBI transcript: **NM_000546.6**
- NCBI protein: **NP_000537.3**
- UniProt: **P04637 (P53_HUMAN)**

## Repository structure

```text
TP53-DNA-Sequence-Analyzer/
├── Code/
│   ├── Linux/
│   │   ├── file_inventory.txt
│   │   ├── linux_commands_history.txt
│   │   └── tp53_analysis.sh
│   ├── Python/
│   │   └── tp53_analysis.py
│   └── R/
│       └── DNA_analysis.R
├── Data/
│   ├── TP53.fa
│   ├── TP53_CDS.fa
│   ├── TP53_protein.fa
│   ├── TP53_uniprot_reference.fa
│   └── comparative TP53 protein sequences
├── Documentation/
├── Plots/
│   ├── Linux/
│   ├── Python/
│   └── R/
├── Results/
│   ├── Linux/
│   ├── Python/
│   └── R/
├── requirements.txt
└── R_packages.txt
```

## Python

### Requirements

```bash
pip install -r requirements.txt
```

### Run

From the repository root:

```bash
python Code/Python/tp53_analysis.py
```

The script locates the input files from the repository structure, so it does **not** depend on a missing `sequence.fasta` in the current working directory.

### Main analyses

- FASTA parsing and validation
- Nucleotide counts and percentages
- GC/AT content and GC/AT ratio
- GC and AT skew
- Frame-aware ORF detection in all three forward reading frames
- Translation of the annotated CDS
- Exact CDS/protein reference validation
- Protein amino-acid composition
- Molecular weight, theoretical pI, instability index, aromaticity and GRAVY using Biopython
- Static nucleotide and GC/AT plots

## Linux

`linux_commands_history.txt` is retained as an **archive of the original command history**.

For a clean reproducible workflow, use:

```bash
bash Code/Linux/tp53_analysis.sh
```

The script performs deterministic local sequence validation and basic composition/property analysis. Its outputs are written to `Results/Linux/Reproduced/`, so the historical Linux results under `Results/Linux/` are not overwritten. It also writes `run_manifest.txt` containing software versions and SHA256 checksums of the core input files and reproduced outputs.

BLASTP is database- and network-dependent. The repository retains the previously recorded BLASTP results under `Results/Linux/`; these are recorded external results, not a fresh query produced by the reproducibility script.

## R / RStudio

Install the packages listed in `R_packages.txt`.

Run from the repository root:

```bash
Rscript Code/R/DNA_analysis.R
```

or open `Code/R/DNA_analysis.R` in RStudio and run it with the repository as the working project.

The R script automatically searches upward from the working directory for the repository root, then uses:

```text
Data/
Results/R/
Plots/R/
```

instead of the old non-existent `R Bio/` paths.

### Main analyses

- CDS nucleotide composition
- GC/AT content
- Codon counting and percentages
- GC content at codon positions
- Translation of the CDS
- Amino-acid composition
- Protein molecular weight and theoretical pI using the `Peptides` package
- ORF analysis in all three frames
- Comparison of the computational ORF with the annotated TP53 CDS
- Static visualizations

There is intentionally no `p3` plot object. The final R figure uses `p1`, `p2`, `p4`, `p5`, and `p6`.

## Validation results

The repository data contain:

- TP53 CDS: **1182 nt**
- TP53 protein: **393 aa** after removing the terminal stop symbol
- CDS matches transcript coordinates **143–1324**
- Translated CDS matches the UniProt P04637 reference
- Protein reference file matches UniProt P04637

The recorded BLASTP result shows a **100% identity, 393-aa alignment** to P53_HUMAN. Because BLAST is an external database query, this result is kept as recorded evidence rather than silently presented as a fresh query.

## Important interpretation notes

- A sequence-pattern match such as an N-X-S/T motif is **not proof that glycosylation occurs**.
- S/T-P motifs are potential phosphorylation-related sequence patterns, not experimental phosphorylation evidence.
- “Mutation hotspot positions” in the Linux results are positions being mapped/analyzed; they are not derived from a mutation-frequency dataset in this repository.
- ORF detection is frame-aware: the stop codon must occur in the same reading frame as its start codon.
- The annotated CDS is already known from the reference record, so the CDS three-frame analysis is best described as **reading-frame/ORF analysis within the annotated CDS**, not independent gene annotation.

## Outputs

Python writes:

```text
Results/Python/TP53_results.txt
Plots/Python/TP53_nucleotide_composition.png
Plots/Python/TP53_gc_at_content.png
```

R writes:

```text
Results/R/TP53_R_results.txt
Plots/R/p1.png
Plots/R/p2.png
Plots/R/p4.png
Plots/R/p5.png
Plots/R/p6.png
Plots/R/final_plot.png
```

The existing Linux results and visualizations are retained as historical project outputs. The cleaned Linux script writes its regenerated deterministic checks to `Results/Linux/Reproduced/`.

## Reproducibility

The cleaned Python, Linux and R scripts resolve the repository structure independently of the caller's working directory and create their output directories automatically. The original Linux command history is preserved separately for transparency. The Linux entry point writes deterministic checks to `Results/Linux/Reproduced/` and never overwrites the historical Linux result archive. A detailed workflow is documented in `Documentation/REPRODUCIBILITY.md`.

## Author

**Priyanshu Kumar**

Independent BSc Biotechnology bioinformatics project.
