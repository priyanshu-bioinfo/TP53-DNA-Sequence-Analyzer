#!/usr/bin/env bash
set -euo pipefail

# TP53 Linux reproducibility workflow
# -----------------------------------
# Run from any directory with:
#   bash /path/to/TP53-DNA-Sequence-Analyzer/Code/Linux/tp53_analysis.sh
#
# The workflow regenerates only deterministic local outputs under:
#   Results/Linux/Reproduced/
# It deliberately does NOT overwrite the historical Linux results in
# Results/Linux/. Historical BLASTP output is external/database-dependent
# and remains preserved as recorded evidence.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA="$ROOT/Data"
OUT="$ROOT/Results/Linux/Reproduced"

mkdir -p "$OUT"

require_cmd() {
    command -v "$1" >/dev/null 2>&1 || {
        echo "ERROR: required command not found: $1" >&2
        exit 1
    }
}

require_cmd python3

# Pass the absolute repository root explicitly; do not depend on the caller's cwd.
export TP53_PROJECT_ROOT="$ROOT"

python3 - <<'PY'
from pathlib import Path
import hashlib
import os
import platform
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

root = Path(os.environ["TP53_PROJECT_ROOT"])
data = root / "Data"
out = root / "Results" / "Linux" / "Reproduced"

required = [
    "TP53.fa",
    "TP53_CDS.fa",
    "TP53_protein.fa",
    "TP53_uniprot_reference.fa",
]

for name in required:
    path = data / name
    if not path.is_file():
        raise SystemExit(f"ERROR: missing required input: {path}")


def read_one(name):
    records = list(SeqIO.parse(data / name, "fasta"))
    if len(records) != 1:
        raise SystemExit(f"ERROR: {name} must contain exactly one FASTA record.")
    return records[0]

transcript = read_one("TP53.fa")
cds = read_one("TP53_CDS.fa")
protein = read_one("TP53_protein.fa")
uniprot = read_one("TP53_uniprot_reference.fa")

transcript_seq = str(transcript.seq).upper()
cds_seq = str(cds.seq).upper()
protein_seq = str(protein.seq).upper().rstrip("*")
uniprot_seq = str(uniprot.seq).upper().rstrip("*")

for label, seq in [("TP53 transcript", transcript_seq), ("TP53 CDS", cds_seq)]:
    if not seq or any(base not in "ACGT" for base in seq):
        raise SystemExit(f"ERROR: {label} contains invalid DNA characters.")

if len(transcript_seq) != 2512:
    raise SystemExit(f"ERROR: unexpected transcript length: {len(transcript_seq)}")
if len(cds_seq) != 1182:
    raise SystemExit(f"ERROR: unexpected CDS length: {len(cds_seq)}")
if len(protein_seq) != 393 or len(uniprot_seq) != 393:
    raise SystemExit("ERROR: unexpected protein/reference length.")

# ------------------------------------------------------------------
# 1. Nucleotide composition of the complete transcript
# ------------------------------------------------------------------
counts = {base: transcript_seq.count(base) for base in "ACGT"}
length = len(transcript_seq)
gc = (counts["G"] + counts["C"]) / length * 100
at = (counts["A"] + counts["T"]) / length * 100

(out / "TP53_nucleotide_composition.txt").write_text(
    "TP53 NUCLEOTIDE COMPOSITION\n"
    "===========================\n\n"
    + "\n".join(f"{b}: {counts[b]} ({counts[b]/length*100:.2f}%)" for b in "ACGT")
    + f"\nTotal: {length} nt\nGC content: {gc:.2f}%\nAT content: {at:.2f}%\n",
    encoding="utf-8",
)

# ------------------------------------------------------------------
# 2. Protein physicochemical properties (Biopython)
# ------------------------------------------------------------------
analysis = ProteinAnalysis(uniprot_seq)
(out / "TP53_physicochemical.txt").write_text(
    "TP53 PHYSICOCHEMICAL PROPERTIES\n"
    "================================\n\n"
    f"Sequence length: {len(uniprot_seq)} aa\n"
    f"Molecular weight: {analysis.molecular_weight():.2f} Da\n"
    f"Theoretical pI (Biopython): {analysis.isoelectric_point():.2f}\n"
    f"Instability index: {analysis.instability_index():.2f}\n"
    f"Aromaticity: {analysis.aromaticity():.4f}\n"
    f"GRAVY: {analysis.gravy():.4f}\n",
    encoding="utf-8",
)

# ------------------------------------------------------------------
# 3. Direct reference validation
# ------------------------------------------------------------------
annotated_cds = transcript_seq[142:1324]  # 1-based coordinates 143-1324
translated = str(Seq(cds_seq).translate()).rstrip("*")
protein_match = protein_seq == uniprot_seq
cds_match = cds_seq == annotated_cds
translation_match = translated == uniprot_seq

lines = [
    "TP53 REFERENCE VALIDATION",
    "=========================",
    "",
    "PROTEIN vs UniProt P04637",
    "------------------------",
    f"Our protein length: {len(protein_seq)} aa",
    f"Reference length: {len(uniprot_seq)} aa",
    f"Exact sequence match: {protein_match}",
    f"Differences: {sum(a != b for a, b in zip(protein_seq, uniprot_seq)) + abs(len(protein_seq)-len(uniprot_seq))}",
    "",
    "CDS vs NCBI NM_000546.6",
    "----------------------",
    f"Our CDS length: {len(cds_seq)} nt",
    "Transcript-annotated CDS coordinates: 143-1324",
    f"Exact sequence match: {cds_match}",
    f"Differences: {sum(a != b for a, b in zip(cds_seq, annotated_cds)) + abs(len(cds_seq)-len(annotated_cds))}",
    "",
    "CDS TRANSLATION vs UniProt P04637",
    "----------------------------------",
    f"Translated protein length: {len(translated)} aa",
    f"Exact sequence match: {translation_match}",
    "",
    "NOTE",
    "----",
    "The BLASTP results in Results/Linux/ are recorded external results.",
    "They are not overwritten by this local reproducibility workflow.",
]
(out / "TP53_reference_validation.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")

# ------------------------------------------------------------------
# 4. Run manifest + checksums for reproducibility/auditability
# ------------------------------------------------------------------
inputs = [data / n for n in required]
outputs = [out / n for n in [
    "TP53_nucleotide_composition.txt",
    "TP53_physicochemical.txt",
    "TP53_reference_validation.txt",
]]

manifest = [
    "TP53 LINUX REPRODUCIBILITY MANIFEST",
    "====================================",
    f"Project root: {root}",
    f"Python: {sys.version.split()[0]}",
    f"Biopython: {__import__('Bio').__version__}",
    f"Platform: {platform.platform()}",
    "",
    "INPUT CHECKSUMS (SHA256)",
    "-------------------------",
]
for path in inputs:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest.append(f"{digest}  {path.relative_to(root)}")

manifest += ["", "REPRODUCED OUTPUT CHECKSUMS (SHA256)", "------------------------------------"]
for path in outputs:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    manifest.append(f"{digest}  {path.relative_to(root)}")

manifest += [
    "",
    "External database results (BLASTP) are intentionally not regenerated here.",
    "The historical BLAST outputs remain in Results/Linux/ and are treated as recorded evidence.",
]
(out / "run_manifest.txt").write_text("\n".join(manifest) + "\n", encoding="utf-8")

print("TP53 Linux reproducibility workflow completed successfully.")
print(f"Project root: {root}")
print(f"Reproduced outputs: {out}")
print("Historical Results/Linux files were not overwritten.")
PY
