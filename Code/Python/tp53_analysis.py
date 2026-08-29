with open("sequence.fasta","r") as file:
    lines = file.readlines()
    sequence = ""
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip()
print("Total DNA bases:",len(sequence))
a_count = sequence.count("A")
t_count = sequence.count("T")
g_count = sequence.count("G")
c_count = sequence.count("C")

print("A:", a_count)
print("T:", t_count)
print("G:", g_count)
print("C:", c_count)
gc_count = g_count + c_count
gc_percentage = (gc_count / len(sequence)) * 100

print("GC content:", gc_percentage, "%")
at_percentage = ((a_count + t_count) / len(sequence)) * 100

print("AT content:", at_percentage, "%")
gc_at_ratio = gc_count / (a_count + t_count)

print("GC/AT ratio:", gc_at_ratio)
atg_count = sequence.count("ATG")

print("ATG count:", atg_count)
taa_count = sequence.count("TAA")
tag_count = sequence.count("TAG")
tga_count = sequence.count("TGA")

print("TAA count:", taa_count)
print("TAG count:", tag_count)
print("TGA count:", tga_count)
codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]

print("First 10 codons:", codons[:10])
coding_regions = []

start = sequence.find("ATG")

while start != -1:
    stop_positions = []

    for stop in ["TAA", "TAG", "TGA"]:
        position = sequence.find(stop, start + 3)
        if position != -1:
            stop_positions.append(position)

    if stop_positions:
        stop = min(stop_positions)
        coding_region = sequence[start:stop + 3]
        coding_regions.append(coding_region)

    start = sequence.find("ATG", start + 1)

print("Number of coding regions:", len(coding_regions))
for i, region in enumerate(coding_regions):
    print("Coding region", i+1, "length:", len(region))
    longest_region = max(coding_regions, key=len)

print("Longest coding region length:", len(longest_region))
print("Longest coding region:", longest_region)
codon_table = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"
}

protein = ""

for i in range(0, len(longest_region) - 2, 3):
    codon = longest_region[i:i+3]
    protein += codon_table.get(codon, "X")

print("Protein sequence:")
print(protein)
# Protein analysis

protein_length = len(protein)

print("Protein length:", protein_length)

# Remove stop symbol for amino acid analysis
protein_without_stop = protein.replace("*", "")

print("Protein length without stop:", len(protein_without_stop))

# Count some important amino acids
print("Alanine (A):", protein_without_stop.count("A"))
print("Glycine (G):", protein_without_stop.count("G"))
print("Proline (P):", protein_without_stop.count("P"))
print("Leucine (L):", protein_without_stop.count("L"))
# Amino acid composition
amino_acids = "ACDEFGHIKLMNPQRSTVWY"

print("Amino acid composition:")

for aa in amino_acids:
    count = protein_without_stop.count(aa)
    percentage = (count / len(protein_without_stop)) * 100
    print(aa, ":", count, "(", percentage, "%)")
   # Molecular weight of protein

amino_acid_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.10,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.20,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "Y": 181.19
}

molecular_weight = 0

for aa in protein_without_stop:
    molecular_weight += amino_acid_weights[aa]

print("Approximate molecular weight:", molecular_weight, "Da") 
# Basic protein composition

hydrophobic = "AILMFWV"
positive = "KR"
negative = "DE"
polar = "NQSTY"

hydrophobic_count = sum(protein_without_stop.count(aa) for aa in hydrophobic)
positive_count = sum(protein_without_stop.count(aa) for aa in positive)
negative_count = sum(protein_without_stop.count(aa) for aa in negative)
polar_count = sum(protein_without_stop.count(aa) for aa in polar)

print("Hydrophobic amino acids:", hydrophobic_count)
print("Positively charged amino acids:", positive_count)
print("Negatively charged amino acids:", negative_count)
print("Polar amino acids:", polar_count)
# Net charge estimate

net_charge = positive_count - negative_count

print("Estimated net charge:", net_charge)
# Approximate isoelectric point (pI)

pI = 7.0

if positive_count > negative_count:
    pI += 0.5
elif positive_count < negative_count:
    pI -= 0.5

print("Approximate pI:", pI)
# Hydrophobic amino acid percentage

hydrophobic_percentage = (hydrophobic_count / len(protein_without_stop)) * 100

print("Hydrophobic percentage:", hydrophobic_percentage, "%")
# Find hydrophobic regions

window_size = 5

for i in range(len(protein_without_stop) - window_size + 1):
    window = protein_without_stop[i:i + window_size]

    hydrophobic_in_window = sum(
        window.count(aa) for aa in hydrophobic
    )

    if hydrophobic_in_window >= 4:
        print("Hydrophobic region:", i + 1, "-", i + window_size,
              window, "(", hydrophobic_in_window, "/5 )")
        # Nucleotide percentages

total = len(sequence)

a_percentage = (a_count / total) * 100
t_percentage = (t_count / total) * 100
g_percentage = (g_count / total) * 100
c_percentage = (c_count / total) * 100

print("\nNucleotide percentages:")
print("A:", round(a_percentage, 2), "%")
print("T:", round(t_percentage, 2), "%")
print("G:", round(g_percentage, 2), "%")
print("C:", round(c_percentage, 2), "%")
# Reverse complement

complement = str.maketrans("ATGC", "TACG")

reverse_complement = sequence.translate(complement)[::-1]

print("\nReverse complement:")
print(reverse_complement[:100])
# GC skew

gc_skew = (g_count - c_count) / (g_count + c_count)

print("\nGC skew:", round(gc_skew, 4))

# AT skew

at_skew = (a_count - t_count) / (a_count + t_count)

print("AT skew:", round(at_skew, 4))

import matplotlib.pyplot as plt

bases = ["A", "T", "G", "C"]
counts = [a_count, t_count, g_count, c_count]

plt.figure(figsize=(7, 5))
plt.bar(bases, counts)

plt.title("TP53 Nucleotide Composition")
plt.xlabel("Nucleotide")
plt.ylabel("Number of Bases")

plt.tight_layout()
plt.show()

composition = ["GC", "AT"]
values = [gc_count, a_count + t_count]

plt.figure(figsize=(7, 5))
plt.bar(composition, values)

plt.title("TP53 GC and AT Composition")
plt.xlabel("Composition")
plt.ylabel("Number of Bases")

plt.tight_layout()
plt.show()

# Save sequence analysis results

results = f"""TP53 DNA Sequence Analysis

Accession: NM_000546.6
Sequence length: {total} bp

Nucleotide counts:
A: {a_count}
T: {t_count}
G: {g_count}
C: {c_count}

Nucleotide percentages:
A: {a_percentage:.2f}%
T: {t_percentage:.2f}%
G: {g_percentage:.2f}%
C: {c_percentage:.2f}%

GC content: {gc_percentage:.2f}%
AT content: {at_percentage:.2f}%
GC/AT ratio: {gc_at_ratio:.4f}

GC skew: {gc_skew:.4f}
AT skew: {at_skew:.4f}
"""

with open("TP53_results.txt", "w") as file:
    file.write(results)

print("\nResults saved to TP53_results.txt")

from Bio import SeqIO

record = SeqIO.read("sequence.fasta", "fasta")

print("\nBiopython analysis")
print("Sequence ID:", record.id)
print("Sequence length:", len(record.seq))

# Biopython analysis

from Bio.Seq import Seq

dna = Seq(sequence)

print("\n--- Biopython Analysis ---")

print("DNA length:", len(dna))

print("Reverse complement:")
print(dna.reverse_complement())

print("Transcription:")
rna = dna.transcribe()
print(rna)

print("Translation:")

# Trim DNA sequence so its length is a multiple of 3
trimmed_dna = dna[:len(dna) - (len(dna) % 3)]

protein_biopython = trimmed_dna.translate(to_stop=False)

print(protein_biopython)
print("\nProtein length:")
print(len(protein_biopython), "amino acids")

print("\nNumber of stop codons:")
print(protein_biopython.count("*"))