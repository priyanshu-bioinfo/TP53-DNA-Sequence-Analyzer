#!/usr/bin/env python3
"""
TP53 DNA Sequence Analyzer
---------------------------
Reproducible Python component of the Human TP53 sequence-analysis project.

The script:
1. Reads the full NM_000546.6 transcript from Data/TP53.fa.
2. Validates the annotated CDS against Data/TP53_CDS.fa.
3. Calculates nucleotide composition, GC/AT content and skews.
4. Finds ORFs in all three forward reading frames using frame-aware codons.
5. Translates the validated CDS and checks it against the UniProt reference.
6. Calculates protein composition and physicochemical properties with Biopython.
7. Generates and saves two Python plots.
8. Writes Results/Python/TP53_results.txt.

Run from any working directory:
    python Code/Python/tp53_analysis.py
"""

from pathlib import Path
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "Data"
RESULTS_DIR = ROOT / "Results" / "Python"
PLOTS_DIR = ROOT / "Plots" / "Python"

TRANSCRIPT_FILE = DATA_DIR / "TP53.fa"
CDS_FILE = DATA_DIR / "TP53_CDS.fa"
PROTEIN_FILE = DATA_DIR / "TP53_protein.fa"
UNIPROT_FILE = DATA_DIR / "TP53_uniprot_reference.fa"


def read_single_fasta(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")
    return SeqIO.read(path, "fasta")


def find_orfs(sequence: str):
    """Find the first in-frame stop after every ATG in each forward frame."""
    sequence = sequence.upper()
    stop_codons = {"TAA", "TAG", "TGA"}
    records = []

    for frame in range(3):
        frame_seq = sequence[frame:]
        usable_length = len(frame_seq) - (len(frame_seq) % 3)
        frame_seq = frame_seq[:usable_length]
        codons = [frame_seq[i:i + 3] for i in range(0, usable_length, 3)]

        for start_idx, codon in enumerate(codons):
            if codon != "ATG":
                continue

            for stop_idx in range(start_idx + 1, len(codons)):
                if codons[stop_idx] in stop_codons:
                    start_nt = frame + start_idx * 3 + 1
                    stop_nt = frame + (stop_idx + 1) * 3
                    nt_length = stop_nt - start_nt + 1
                    aa_length = stop_idx - start_idx
                    records.append(
                        {
                            "frame": frame + 1,
                            "start_codon_index": start_idx + 1,
                            "stop_codon_index": stop_idx + 1,
                            "start_nt": start_nt,
                            "stop_nt": stop_nt,
                            "stop_codon": codons[stop_idx],
                            "length_nt": nt_length,
                            "protein_length_aa": aa_length,
                            "sequence": frame_seq[start_idx * 3:(stop_idx + 1) * 3],
                        }
                    )
                    break

    return records


def format_orf(record):
    return (
        f"Frame {record['frame']}: {record['start_nt']}-{record['stop_nt']} nt; "
        f"{record['length_nt']} nt / {record['protein_length_aa']} aa; "
        f"stop={record['stop_codon']}"
    )


def main():
    parser = argparse.ArgumentParser(description="Analyze human TP53 sequence data.")
    parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    transcript_record = read_single_fasta(TRANSCRIPT_FILE)
    cds_record = read_single_fasta(CDS_FILE)
    protein_record = read_single_fasta(PROTEIN_FILE)
    uniprot_record = read_single_fasta(UNIPROT_FILE)

    transcript = str(transcript_record.seq).upper()
    cds = str(cds_record.seq).upper()
    protein_with_stop = str(protein_record.seq).upper()
    protein = protein_with_stop.rstrip("*")
    uniprot = str(uniprot_record.seq).upper().rstrip("*")

    if set(transcript) - set("ACGT"):
        raise ValueError("Transcript contains characters outside A/C/G/T.")
    if set(cds) - set("ACGT"):
        raise ValueError("CDS contains characters outside A/C/G/T.")

    # Nucleotide analysis
    counts = {base: transcript.count(base) for base in "ATGC"}
    total = len(transcript)
    gc_count = counts["G"] + counts["C"]
    at_count = counts["A"] + counts["T"]
    gc_pct = gc_count / total * 100
    at_pct = at_count / total * 100
    gc_at_ratio = gc_count / at_count
    gc_skew = (counts["G"] - counts["C"]) / gc_count
    at_skew = (counts["A"] - counts["T"]) / at_count

    # CDS validation
    expected_cds = transcript[142:1324]  # NCBI NM_000546.6 coordinates 143-1324
    cds_matches_transcript = cds == expected_cds

    translated_cds = str(Seq(cds).translate()).rstrip("*")
    protein_matches_uniprot = translated_cds == uniprot
    protein_file_matches_uniprot = protein == uniprot

    # ORF analysis on the full transcript
    orfs = find_orfs(transcript)
    longest_by_frame = {}
    for frame in (1, 2, 3):
        frame_orfs = [o for o in orfs if o["frame"] == frame]
        if frame_orfs:
            longest_by_frame[frame] = max(frame_orfs, key=lambda x: x["length_nt"])

    if not longest_by_frame:
        raise RuntimeError("No ORFs were found in the transcript.")

    longest_orf = max(longest_by_frame.values(), key=lambda x: x["length_nt"])

    # Protein properties on the canonical 393-aa protein.
    analysis = ProteinAnalysis(uniprot)
    aa_counts = analysis.count_amino_acids()
    aa_percent = analysis.amino_acids_percent

    # Results file
    out = []
    out.append("TP53 DNA SEQUENCE ANALYSIS — PYTHON")
    out.append("=" * 40)
    out.append("")
    out.append(f"Transcript: {transcript_record.id}")
    out.append(f"Transcript length: {total} nt")
    out.append("")
    out.append("FULL TRANSCRIPT NUCLEOTIDE COMPOSITION")
    out.append("----------------------------------------")
    for base in "ATGC":
        out.append(
            f"{base}: {counts[base]} ({counts[base] / total * 100:.2f}%)"
        )
    out.append(f"GC content: {gc_pct:.2f}%")
    out.append(f"AT content: {at_pct:.2f}%")
    out.append(f"GC/AT ratio: {gc_at_ratio:.4f}")
    out.append(f"GC skew: {gc_skew:.4f}")
    out.append(f"AT skew: {at_skew:.4f}")
    out.append("")
    out.append("CDS / PROTEIN VALIDATION")
    out.append("------------------------")
    out.append(f"CDS length: {len(cds)} nt")
    out.append(f"CDS matches transcript positions 143-1324: {cds_matches_transcript}")
    out.append(f"Translated CDS length: {len(translated_cds)} aa")
    out.append(f"Protein file length: {len(protein)} aa")
    out.append(f"UniProt reference length: {len(uniprot)} aa")
    out.append(f"Translated CDS matches UniProt P04637: {protein_matches_uniprot}")
    out.append(f"Protein file matches UniProt P04637: {protein_file_matches_uniprot}")
    out.append("")
    out.append("FULL-TRANSCRIPT ORF ANALYSIS")
    out.append("----------------------------")
    for frame in (1, 2, 3):
        if frame in longest_by_frame:
            out.append(format_orf(longest_by_frame[frame]))
        else:
            out.append(f"Frame {frame}: no complete ORF found")
    out.append(f"Longest ORF overall: {format_orf(longest_orf)}")
    out.append("")
    out.append("CANONICAL TP53 PROTEIN PROPERTIES")
    out.append("----------------------------------")
    out.append(f"Protein length: {len(uniprot)} aa")
    out.append(f"Molecular weight (Biopython): {analysis.molecular_weight():.2f} Da")
    out.append(f"Theoretical pI (Biopython): {analysis.isoelectric_point():.2f}")
    out.append(f"Instability index: {analysis.instability_index():.2f}")
    out.append(f"Aromaticity: {analysis.aromaticity():.4f}")
    out.append(f"GRAVY: {analysis.gravy():.4f}")
    out.append("")
    out.append("AMINO ACID COMPOSITION")
    out.append("----------------------")
    for aa in "ACDEFGHIKLMNPQRSTVWY":
        out.append(f"{aa}: {aa_counts[aa]} ({aa_percent[aa]:.2f}%)")
    out.append("")
    out.append("NOTES")
    out.append("-----")
    out.append("ORFs are detected only when the stop codon is in the same reading frame as the ATG.")
    out.append("Protein pI and physicochemical properties are calculated with Biopython ProteinAnalysis.")
    out.append("A motif match is a sequence-pattern observation and is not evidence of biological modification.")
    (RESULTS_DIR / "TP53_results.txt").write_text("\n".join(out) + "\n", encoding="utf-8")

    # Plot 1: nucleotide composition
    plt.figure(figsize=(7, 5))
    plt.bar(list("ATGC"), [counts[b] for b in "ATGC"])
    plt.title("Human TP53 Nucleotide Composition")
    plt.xlabel("Nucleotide")
    plt.ylabel("Number of Bases")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "TP53_nucleotide_composition.png", dpi=300)
    plt.close()

    # Plot 2: GC vs AT
    plt.figure(figsize=(7, 5))
    plt.bar(["GC", "AT"], [gc_pct, at_pct])
    plt.title("Human TP53 GC and AT Content")
    plt.xlabel("Composition")
    plt.ylabel("Percentage (%)")
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "TP53_gc_at_content.png", dpi=300)
    plt.close()

    print(f"Analysis complete. Results: {RESULTS_DIR / 'TP53_results.txt'}")
    print(f"Plots: {PLOTS_DIR}")


if __name__ == "__main__":
    main()
