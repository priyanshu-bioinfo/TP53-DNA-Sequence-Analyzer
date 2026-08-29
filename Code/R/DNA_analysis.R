# Save all plots
dir.create("R Bio/Plots", recursive = TRUE, showWarnings = FALSE)

ggsave("R Bio/Plots/p1.png", p1, width = 10, height = 7)
ggsave("R Bio/Plots/p2.png", p2, width = 10, height = 7)
ggsave("R Bio/Plots/p4.png", p4, width = 10, height = 7)
ggsave("R Bio/Plots/p5.png", p5, width = 10, height = 7)
ggsave("R Bio/Plots/p6.png", p6, width = 10, height = 7)

ggsave("R Bio/Plots/final_plot.png", final_plot, width = 12, height = 10)



# DNA Sequence Analyzer
# Target gene: TP53

library(Biostrings)

tp53_dna <- readDNAStringSet("R Bio/Data/sequence.fasta")
tp53_dna
tp53_length <- width(tp53_dna)
tp53_length
tp53_counts <-alphabetFrequency(tp53_dna)
tp53_counts

# GC content
gc_count <- letterFrequency(tp53_dna, letters = c("G", "C"))

gc_content <- sum(gc_count) / tp53_length * 100

gc_content

tp53_percent <- (tp53_counts[, c("A", "C", "G", "T")] / tp53_length) * 100
tp53_percent

# Calculate AT content

at_count <- letterFrequency(tp53_dna, letters = c("A", "T"))

at_content <- sum(at_count) / tp53_length * 100

at_content

# Plot TP53 nucleotide composition

barplot(
  tp53_counts[1,c("A", "C", "G", "T")],
  main = "TP53 Nucleotide Composition",
  xlab = "Nucleotide",
  ylab = "Count"
)

# Plot TP53 nucleotide percentages

barplot(
  tp53_percent,
  main = "TP53 Nucleotide Percentage",
  xlab = "Nucleotide",
  ylab = "Percentage (%)"
)

# Plot GC and AT content

gc_at <- c(
  GC = gc_content,
  AT = at_content
)

barplot(
  gc_at,
  main = "TP53 GC and AT Content",
  xlab = "Base Composition",
  ylab = "Percentage (%)"
)

# Convert TP53 DNA sequence to RNA sequence

tp53_rna <- RNAStringSet(tp53_dna)

tp53_rna

class(tp53_rna)

width(tp53_rna)

alphabetFrequency(tp53_rna, baseOnly = TRUE)

# Translate TP53 DNA sequence into protein
tp53_protein <- translate(tp53_dna)

tp53_protein

tp53_cds <- subseq(tp53_dna, start = 143, end = 1324)
width(tp53_cds)
tp53_protein <- translate(tp53_cds)
tp53_protein
width(tp53_protein)

# Translate the TP53 coding sequence into protein
tp53_protein_cds <- translate(tp53_cds)

# Check the protein length
width(tp53_protein_cds)

# Display the translated protein
tp53_protein_cds

# Remove the final stop codon from the translated TP53 protein
tp53_protein_final <- subseq(
  tp53_protein_cds,
  start = 1,
  end = 393
)
width(tp53_protein_final)
tp53_protein_final

# Count amino acids in the TP53 protein
tp53_aa_counts <- alphabetFrequency(
  tp53_protein_final,
  baseOnly = TRUE
)

tp53_aa_counts
# Calculate the percentage of each amino acid
tp53_aa_percent <- (tp53_aa_counts / width(tp53_protein_final)) * 100

tp53_aa_percent

# Plot TP53 amino-acid composition

barplot(
  as.vector(tp53_aa_percent),
  names.arg = colnames(tp53_aa_percent),
  main = "TP53 Amino Acid Composition",
  xlab = "Amino Acid",
  ylab = "Percentage (%)",
  las = 2
)

# Calculate GC content of the TP53 coding sequence

cds_counts <- letterFrequency(
  tp53_cds,
  letters = c("G", "C")
)

cds_gc_content <- sum(cds_counts) / width(tp53_cds) * 100

cds_gc_content

# Calculate AT content of the TP53 coding sequence

cds_at_counts <- letterFrequency(
  tp53_cds,
  letters = c("A", "T")
)

cds_at_content <- sum(cds_at_counts) / width(tp53_cds) * 100

cds_at_content

# Plot CDS GC and AT content

cds_gc_at <- c(
  GC = cds_gc_content,
  AT = cds_at_content
)

barplot(
  cds_gc_at,
  main = "TP53 CDS GC and AT Content",
  xlab = "Base Composition",
  ylab = "Percentage (%)"
)

# Analyze TP53 codon usage

tp53_codons <- oligonucleotideFrequency(
  tp53_cds,
  width = 3
)
tp53_codons

# Find the 10 most frequent TP53 codons

top_10_codons <- sort(
  tp53_codons,
  decreasing = TRUE
)[1:10]
top_10_codons

# Plot the 10 most frequent TP53 codons

barplot(
  top_10_codons,
  main = "Top 10 Most Frequent TP53 Codons",
  xlab = "Codon",
  ylab = "Frequency",
  las = 2
)

# Reading Frame 1
frame1_dna <- subseq(
  tp53_dna,
  start = 1,
  end = 2511
)

frame1_protein <- translate(frame1_dna)

# Check the length of Frame 1
width(frame1_dna)
width(frame1_protein)

# Display Frame 1 protein
frame1_protein

# Reading Frame 2
frame2_dna <- subseq(
  tp53_dna,
  start = 2,
  end = 2512
)
frame2_protein <- translate(frame2_dna)

# Check the length of Frame 2
width(frame2_dna)
width(frame2_protein)

# Display Frame 2 protein
frame2_protein

# Reading Frame 3
frame3_dna <- subseq(
  tp53_dna,
  start = 3,
  end = 2510
)

frame3_protein <- translate(frame3_dna)

# Check the length of Frame 3
width(frame3_dna)
width(frame3_protein)

# Display Frame 3 protein
frame3_protein

# ORF analysis - Frame 1
frame1_sequence <- as.character(frame1_dna[[1]])
frame1_codons <- substring(
  frame1_sequence,
  seq(1, nchar(frame1_sequence) - 2, by = 3),
  seq(3, nchar(frame1_sequence), by = 3)
)

# Find start codons (ATG)
frame1_start_positions <- which(frame1_codons == "ATG")

# Find stop codons (TAA, TAG, TGA)
frame1_stop_positions <- which(
  frame1_codons %in% c("TAA", "TAG", "TGA")
)
frame1_start_positions
frame1_stop_positions

# ORF analysis - Frame 2
frame2_sequence <- as.character(frame2_dna[[1]])
frame2_codons <- substring(
  frame2_sequence,
  seq(1, nchar(frame2_sequence) - 2, by = 3),
  seq(3, nchar(frame2_sequence), by = 3)
)

# Find start codons (ATG)
frame2_start_positions <- which(frame2_codons == "ATG")

# Find stop codons (TAA, TAG, TGA)
frame2_stop_positions <- which(
  frame2_codons %in% c("TAA", "TAG", "TGA")
)
frame2_start_positions
frame2_stop_positions

# ORF analysis - Frame 3
frame3_sequence <- as.character(frame3_dna[[1]])
frame3_codons <- substring(
  frame3_sequence,
  seq(1, nchar(frame3_sequence) - 2, by = 3),
  seq(3, nchar(frame3_sequence), by = 3)
)

# Find start codons (ATG)
frame3_start_positions <- which(frame3_codons == "ATG")

# Find stop codons (TAA, TAG, TGA)
frame3_stop_positions <- which(
  frame3_codons %in% c("TAA", "TAG", "TGA")
)
frame3_start_positions
frame3_stop_positions

# Function to identify candidate ORFs in a reading frame
find_orfs <- function(start_positions, stop_positions) {
  
  orfs <- data.frame(
    start_codon = integer(),
    stop_codon = integer(),
    codon_length = integer(),
    nucleotide_length = integer()
  )
  
  for (start in start_positions) {
    
    possible_stops <- stop_positions[stop_positions > start]
    
    if (length(possible_stops) > 0) {
      
      stop <- possible_stops[1]
      
      orfs <- rbind(
        orfs,
        data.frame(
          start_codon = start,
          stop_codon = stop,
          codon_length = stop - start + 1,
          nucleotide_length = (stop - start + 1) * 3
        )
      )
    }
  }
  
  return(orfs)
}

# Find candidate ORFs in all three reading frames
frame1_orfs <- find_orfs(
  frame1_start_positions,
  frame1_stop_positions
)
frame2_orfs <- find_orfs(
  frame2_start_positions,
  frame2_stop_positions
)
frame3_orfs <- find_orfs(
  frame3_start_positions,
  frame3_stop_positions
)

# Display the ORFs found in each frame
frame1_orfs
frame2_orfs
frame3_orfs

# Combine ORFs from all three reading frames
frame1_orfs$frame <- "Frame 1"
frame2_orfs$frame <- "Frame 2"
frame3_orfs$frame <- "Frame 3"

all_orfs <- rbind(
  frame1_orfs,
  frame2_orfs,
  frame3_orfs
)
all_orfs

# Find the longest candidate ORF
longest_orf <- all_orfs[
  which.max(all_orfs$nucleotide_length),
]
longest_orf

# NCBI reference annotation for TP53 NM_000546.6
ncbi_cds_start <- 143
ncbi_cds_end <- 1324

ncbi_cds_length <- ncbi_cds_end - ncbi_cds_start + 1

ncbi_cds_length

# Compare the NCBI CDS with our Frame 2 ORF analysis
ncbi_cds_frame <- "Frame 2"

ncbi_start_codon_position <- ((ncbi_cds_start - 2) / 3) + 1

ncbi_stop_codon_position <- ((ncbi_cds_end - 2) / 3) + 1

ncbi_start_codon_position
ncbi_stop_codon_position

# Verify the NCBI-annotated TP53 CDS
reference_tp53_protein <- translate(tp53_cds)

reference_protein_length <- width(reference_tp53_protein) - 1

reference_protein_length

# Compare computational ORFs with the NCBI-annotated CDS
ncbi_orf_match <- all_orfs[
  all_orfs$frame == "Frame 2" &
    all_orfs$start_codon == ncbi_start_codon_position &
    all_orfs$stop_codon == ncbi_stop_codon_position,
]

ncbi_orf_match

names(all_orfs)
head(all_orfs)

all_orfs[all_orfs$frame == "Frame 2", ]
table(all_orfs$frame)
frame2_orfs <- all_orfs[all_orfs$frame == "Frame 2", ]
frame2_orfs
aggregate(nucleotide_length ~frame, data = all_orfs, max)
all_orfs[all_orfs$nucleotide_length ==
           max(all_orfs$nucleotide_length),]

# Extract the longest computational ORF
longest_orf <- all_orfs[all_orfs$nucleotide_length ==
                          max(all_orfs$nucleotide_length), ]

longest_orf
# Extract the NCBI-annotated CDS
ncbi_cds_seq <- subseq(tp53_dna, start = 143, end = 1324)

width(ncbi_cds_seq)

# Extract the computational ORF nucleotide sequence
orf_start_nt <- ((longest_orf$start_codon - 1) * 3) + 2
orf_end_nt <- (longest_orf$stop_codon * 3) + 1

computational_orf_seq <- subseq(
  tp53_dna,
  start = orf_start_nt,
  end = orf_end_nt
)

width(computational_orf_seq)
# Compare computational ORF with NCBI CDS
computational_orf_seq == ncbi_cds_seq

# Translate the computational ORF
computational_orf_protein <- translate(computational_orf_seq)

# Translate the NCBI CDS
ncbi_cds_protein <- translate(ncbi_cds_seq)

# Check protein lengths
width(computational_orf_protein)
width(ncbi_cds_protein)

# Compare the translated proteins
computational_orf_protein == ncbi_cds_protein
computational_orf_protein

# Nucleotide composition of the NCBI-annotated TP53 CDS
seq_chars <- strsplit(as.character(ncbi_cds_seq), "")[[1]]

nucleotide_counts <- table(seq_chars)

nucleotide_counts

# Calculate percentages
nucleotide_percent <- prop.table(nucleotide_counts) * 100

nucleotide_percent

# GC content
gc_content <- sum(nucleotide_counts[c("G", "C")]) /
  sum(nucleotide_counts) * 100
gc_content

# Nucleotide composition of the complete TP53 DNA sequence
full_seq_chars <- strsplit(as.character(tp53_dna), "")[[1]]
full_nucleotide_counts <- table(full_seq_chars)
full_nucleotide_counts

# Percentages
full_nucleotide_percent <- prop.table(full_nucleotide_counts) * 100
full_nucleotide_percent

# GC content of the complete sequence
full_gc_content <- sum(full_nucleotide_counts[c("G", "C")]) /
  sum(full_nucleotide_counts) * 100
full_gc_content

# Longest ORF from each reading frame
longest_by_frame <- do.call(
  rbind,
  lapply(split(all_orfs, all_orfs$frame), function(x) {
    x[which.max(x$ nucleotide_length), ]
  })
)

longest_by_frame

longest_orf <-longest_by_frame["Frame 2", ]
longest_orf

# Convert ORF codon positions to nucleotide positions
orf_start_nt <- (longest_orf$start_codon - 1) * 3 + 2
orf_end_nt <- (longest_orf$stop_codon - 1) * 3 + 2 + 2

computational_orf_seq <- subseq(
  tp53_dna,
  start = orf_start_nt,
  end = orf_end_nt
)
width(computational_orf_seq)

# Translate the computational ORF
computational_orf_protein <- translate(computational_orf_seq)

# Translate the NCBI CDS
ncbi_cds_protein <- translate(ncbi_cds_seq)

# Check protein lengths
width(computational_orf_protein)
width(ncbi_cds_protein)

# Compare the translated proteins
computational_orf_protein == ncbi_cds_protein

# Nucleotide composition of the NCBI-annotated TP53 CDS
seq_chars <- strsplit(as.character(ncbi_cds_seq), "")[[1]]

nucleotide_counts <- table(seq_chars)

nucleotide_counts

# Calculate percentages
nucleotide_percent <- prop.table(nucleotide_counts) * 100

nucleotide_percent

# GC content
gc_content <- sum(nucleotide_counts[c("G", "C")]) /
  sum(nucleotide_counts) * 100

gc_content

# Nucleotide composition of the complete TP53 DNA sequence
full_seq_chars <- strsplit(as.character(tp53_dna), "")[[1]]

full_nucleotide_counts <- table(full_seq_chars)

full_nucleotide_counts

# Percentages
full_nucleotide_percent <- prop.table(full_nucleotide_counts) * 100

full_nucleotide_percent

# GC content of the complete sequence
full_gc_content <- sum(full_nucleotide_counts[c("G", "C")]) /
  sum(full_nucleotide_counts) * 100

full_gc_content

# Summary of ORFs in the three reading frames
orf_summary <- data.frame(
  Frame = c("Frame 1", "Frame 2", "Frame 3"),
  Number_of_ORFs = c(
    sum(all_orfs$frame == "Frame 1"),
    sum(all_orfs$frame == "Frame 2"),
    sum(all_orfs$frame == "Frame 3")
  ),
  Longest_ORF_nt = c(
    max(all_orfs$nucleotide_length[all_orfs$frame == "Frame 1"]),
    max(all_orfs$nucleotide_length[all_orfs$frame == "Frame 2"]),
    max(all_orfs$nucleotide_length[all_orfs$frame == "Frame 3"])
  )
)
orf_summary

# Plot nucleotide composition of the complete TP53 sequence
barplot(
  full_nucleotide_percent[c("A", "T", "G", "C")],
  main = "TP53 Nucleotide Composition",
  xlab = "Nucleotide",
  ylab = "Percentage (%)"
)

orf_summary

# Plot the longest ORF from each reading frame
barplot(
  orf_summary$Longest_ORF_nt,
  names.arg = orf_summary$Frame,
  main = "Longest ORF in Each Reading Frame",
  xlab = "Reading Frame",
  ylab = "Length (nt)"
)

# Plot the longest ORF from each reading frame
barplot(
  orf_summary$Longest_ORF_nt,
  names.arg = orf_summary$Frame,
  main = "Longest ORF in Each Reading Frame",
  xlab = "Reading Frame",
  ylab = "Length (nt)"
)

# GC content of the computational ORF
orf_chars <- strsplit(as.character(computational_orf_seq), "")[[1]]
orf_counts <- table(orf_chars)
orf_gc_content <- sum(orf_counts[c("G", "C")]) /
  sum(orf_counts) * 100
orf_gc_content

# Compare computational ORF with NCBI CDS
identical(
  as.character(computational_orf_seq),
  as.character(ncbi_cds_seq)
)

# Compare computationally translated protein with NCBI protein
identical(
  as.character(computational_orf_protein),
  as.character(ncbi_cds_protein)
)
# Amino-acid composition of TP53 protein
protein_chars <- strsplit(
  as.character(ncbi_cds_protein),
  ""
)[[1]]
amino_acid_counts <- table(protein_chars)
amino_acid_counts

# Remove the stop symbol from the protein sequence
protein_chars_no_stop <- protein_chars[protein_chars != "*"]

# Calculate amino-acid percentages
amino_acid_percent <- prop.table(
  table(protein_chars_no_stop)
) * 100

amino_acid_percent

# Plot amino-acid composition of TP53 protein
barplot(
  amino_acid_percent,
  main = "TP53 Amino-Acid Composition",
  xlab = "Amino Acid",
  ylab = "Percentage (%)",
  las = 2
)

# Codon composition of the TP53 CDS
cds_sequence <- as.character(ncbi_cds_seq)
codons <- substring(
  cds_sequence,
  seq(1, nchar(cds_sequence) - 2, by = 3),
  seq(3, nchar(cds_sequence), by = 3)
)
codon_counts <- table(codons)
codon_counts

# Remove the terminal stop codon
codon_counts_no_stop <- codon_counts[
  !names(codon_counts) %in% c("TAA", "TAG", "TGA")
]

# Calculate codon percentages
codon_percent <- prop.table(codon_counts_no_stop) * 100
codon_percent

# Plot codon usage of TP53 CDS
barplot(
  codon_percent,
  main = "TP53 Codon Usage",
  xlab = "Codon",
  ylab = "Percentage (%)",
  las = 2,
  cex.names = 0.7
)

# Top 10 most frequent codons in TP53 CDS
top_10_codons <- sort(
  codon_counts_no_stop,
  decreasing = TRUE
)[1:10]

top_10_codons

# Calculate Relative Synonymous Codon Usage (RSCU)
codon_amino_acid <- translate(
  DNAString(
    paste0(
      names(codon_counts_no_stop),
      collapse = ""
    )
  )
)

codon_table <- data.frame(
  Codon = names(codon_counts_no_stop),
  Count = as.numeric(codon_counts_no_stop),
  Amino_Acid = as.character(codon_amino_acid)
)
codon_table$RSCU <- ave(
  codon_table$Count,
  codon_table$Amino_Acid,
  FUN = function(x) x / mean(x)
)
codon_table

# Most preferred codons by RSCU
top_rscu <- codon_table[
  order(codon_table$RSCU, decreasing = TRUE),
]
top_rscu[1:10, ]

# Total number of sense codons
total_sense_codons <- sum(codon_counts_no_stop)
total_sense_codons

# GC content at each codon position
sense_codons <- codons[
  !codons %in% c("TAA", "TAG", "TGA")
]
codon_position_1 <- substring(sense_codons, 1, 1)
codon_position_2 <- substring(sense_codons, 2, 2)
codon_position_3 <- substring(sense_codons, 3, 3)

gc1 <- mean(codon_position_1 %in% c("G", "C")) * 100
gc2 <- mean(codon_position_2 %in% c("G", "C")) * 100
gc3 <- mean(codon_position_3 %in% c("G", "C")) * 100

gc1
gc2
gc3

# Plot GC content at each codon position
gc_by_position <- c(
  GC1 = gc1,
  GC2 = gc2,
  GC3 = gc3
)

barplot(
  gc_by_position,
  main = "TP53 GC Content by Codon Position",
  xlab = "Codon Position",
  ylab = "GC Content (%)"
)

# Protein properties of TP53
tp53_protein <- ncbi_cds_protein

# Remove the stop symbol
tp53_protein_char <- as.character(tp53_protein)

tp53_protein_no_stop <- gsub("\\*", "", tp53_protein_char)

nchar(tp53_protein_no_stop)

# install.packages("Peptides")
library(Peptides)
# Calculate molecular weight
tp53_mw <- mw(tp53_protein_no_stop)
tp53_mw
# Calculate theoretical isoelectric point
tp53_pI <- pI(tp53_protein_no_stop)
tp53_pI

unique(strsplit(as.character(tp53_protein_no_stop),"")[[1]])
tp53_pI
table(strsplit(as.character(tp53_protein_no_stop),"")[[1]])
# Amino acid composition of TP53 protein (%)

aa_counts <- table(strsplit(as.character(tp53_protein_no_stop), "")[[1]])

aa_percent <- (aa_counts / sum(aa_counts)) * 100

aa_percent

# Plot aa composition of TP53 protein
barplot(aa_percent,main = "TP53 Amino Acid Composition",
        xlab="Amino Acid",
        ylab="Percentage(%)"
        )
nchar(as.character(cds_sequence))
# Nucleotide composition of TP53 CDS

cds_chars <- strsplit(as.character(cds_sequence), "")[[1]]

nt_counts <- table(cds_chars)

nt_percent <- (nt_counts / sum(nt_counts)) * 100

nt_counts
nt_percent

# GC content of TP53 CDS

gc_content <- ((nt_counts["G"] + nt_counts["C"]) / sum(nt_counts)) * 100
gc_content

# AT content of TP53 CDS

at_content <- ((nt_counts["A"] + nt_counts["T"]) / sum(nt_counts)) * 100
at_content
# Compare our CDS length with NCBI CDS length
cds_length <- nchar(as.character(cds_sequence))
cds_length
ncbi_cds_length
# GC content of NCBI CDS
ncbi_cds_chars <- strsplit(as.character(ncbi_cds_seq), "")[[1]]
ncbi_nt_counts <- table(ncbi_cds_chars)
ncbi_gc_content <- ((ncbi_nt_counts["G"] + ncbi_nt_counts["C"]) /
                      sum(ncbi_nt_counts)) * 100
ncbi_gc_content

# AT content of NCBI CDS
ncbi_at_content <- ((ncbi_nt_counts["A"] + ncbi_nt_counts["T"]) /
                      sum(ncbi_nt_counts)) * 100

ncbi_at_content

# Protein length comparison
our_protein_length <- nchar(as.character(tp53_protein_no_stop))
ncbi_protein_length <- nchar(as.character(ncbi_cds_protein))
our_protein_length
ncbi_protein_length

# Remove stop symbol from NCBI protein
ncbi_protein_no_stop <- gsub("\\*", "", as.character(ncbi_cds_protein))
ncbi_protein_length <- nchar(ncbi_protein_no_stop)
ncbi_protein_length

# Compare our protein length with NCBI protein length
protein_length <- nchar(as.character(tp53_protein_no_stop))
protein_length
ncbi_protein_length <- nchar(as.character(ncbi_protein_no_stop))
ncbi_protein_length

ls(pattern = "ncbi")

# Compare molecular weight of our protein with NCBI protein
our_mw <- mw(tp53_protein_no_stop)
ncbi_mw <- mw(ncbi_protein_no_stop)
our_mw
ncbi_mw

# Compare theoretical pI of our protein with NCBI protein
our_pI <- pI(tp53_protein_no_stop)
ncbi_pI <- pI(ncbi_protein_no_stop)
our_pI
ncbi_pI

# Nucleotide composition of NCBI CDS
ncbi_nt_chars <- strsplit(as.character(ncbi_cds_seq), "")[[1]]
ncbi_nt_counts <- table(ncbi_nt_chars)
ncbi_nt_percent <- (ncbi_nt_counts / sum(ncbi_nt_counts)) * 100
ncbi_nt_counts
ncbi_nt_percent

# Compare nucleotide composition: our data vs NCBI
nucleotide_comparison <- data.frame(
  Nucleotide = c("A", "C", "G", "T"),
  Our_Data = as.numeric(nt_percent[c("A", "C", "G", "T")]),
  NCBI = as.numeric(ncbi_nt_percent[c("A", "C", "G", "T")])
)

nucleotide_comparison$Difference <-
  nucleotide_comparison$Our_Data - nucleotide_comparison$NCBI

nucleotide_comparison

# Create the three forward reading frames
frame1 <- substring(
  as.character(cds_sequence),
  1,
  nchar(cds_sequence)
)

frame2 <- substring(
  as.character(cds_sequence),
  2,
  nchar(cds_sequence)
)

frame3 <- substring(
  as.character(cds_sequence),
  3,
  nchar(cds_sequence)
)

# Make each frame divisible by 3
frame1 <- substr(frame1, 1, nchar(frame1) - nchar(frame1) %% 3)
frame2 <- substr(frame2, 1, nchar(frame2) - nchar(frame2) %% 3)
frame3 <- substr(frame3, 1, nchar(frame3) - nchar(frame3) %% 3)

# Display frame lengths
nchar(frame1)
nchar(frame2)
nchar(frame3)

# Function to find ORFs in a reading frame
find_orfs <- function(sequence) {
  
  seq <- toupper(sequence)
  
  # Split sequence into codons
  codons <- substring(
    seq,
    seq(1, nchar(seq) - 2, by = 3),
    seq(3, nchar(seq), by = 3)
  )
  
  start_positions <- which(codons == "ATG")
  
  orfs <- data.frame(
    Start_Codon = integer(),
    Stop_Codon = character(),
    Start_Position = integer(),
    Stop_Position = integer(),
    Length_nt = integer(),
    Protein_Length_aa = integer()
  )
  
  for (start in start_positions) {
    
    stops <- which(
      codons %in% c("TAA", "TAG", "TGA") &
        seq_along(codons) > start
    )
    
    if (length(stops) > 0) {
      
      stop <- stops[1]
      
      orfs <- rbind(
        orfs,
        data.frame(
          Start_Codon = start,
          Stop_Codon = codons[stop],
          Start_Position = (start - 1) * 3 + 1,
          Stop_Position = stop * 3,
          Length_nt = (stop - start + 1) * 3,
          Protein_Length_aa = stop - start
        )
      )
    }
  }
  
  return(orfs)
}

# Find ORFs in all three reading frames
orfs_frame1 <- find_orfs(frame1)
orfs_frame2 <- find_orfs(frame2)
orfs_frame3 <- find_orfs(frame3)

# Display the results
orfs_frame1
orfs_frame2
orfs_frame3

# Find the longest ORF in each reading frame
longest_orf_frame1 <- orfs_frame1[
  which.max(orfs_frame1$Length_nt),
]

longest_orf_frame2 <- orfs_frame2[
  which.max(orfs_frame2$Length_nt),
]

longest_orf_frame3 <- orfs_frame3[
  which.max(orfs_frame3$Length_nt),
]

# Display the longest ORF from each frame
longest_orf_frame1
longest_orf_frame2
longest_orf_frame3

# Compare longest ORF length from each reading frame

orf_lengths <- c(
  Frame1 = longest_orf_frame1$Length_nt,
  Frame2 = longest_orf_frame2$Length_nt,
  Frame3 = longest_orf_frame3$Length_nt
)

orf_lengths

# Plot longest ORF length
barplot(
  orf_lengths,
  main = "Longest ORF Length in Each Reading Frame",
  xlab = "Reading Frame",
  ylab = "ORF Length (nt)"
)

# Compare protein lengths of the longest ORFs

orf_protein_lengths <- c(
  Frame1 = longest_orf_frame1$Protein_Length_aa,
  Frame2 = longest_orf_frame2$Protein_Length_aa,
  Frame3 = longest_orf_frame3$Protein_Length_aa
)

orf_protein_lengths

# Plot longest ORF protein lengths
barplot(
  orf_protein_lengths,
  main = "Longest ORF Protein Length",
  xlab = "Reading Frame",
  ylab = "Protein Length (aa)"
)

# Final summary of TP53 DNA sequence analysis
cat("TP53 DNA SEQUENCE ANALYSIS SUMMARY\n")
cat("----------------------------------\n")
cat("CDS length:", cds_length, "nt\n")
cat("Protein length:", protein_length, "aa\n")
cat("GC content:", gc_content, "%\n")
cat("AT content:", at_content, "%\n")
cat("Molecular weight:", our_mw, "Da\n")
cat("Theoretical pI:", our_pI, "\n")
cat("Longest ORF - Frame 1:", longest_orf_frame1$Length_nt, "nt /",
    longest_orf_frame1$Protein_Length_aa, "aa\n")
cat("Longest ORF - Frame 2:", longest_orf_frame2$Length_nt, "nt /",
    longest_orf_frame2$Protein_Length_aa, "aa\n")
cat("Longest ORF - Frame 3:", longest_orf_frame3$Length_nt, "nt /",
    longest_orf_frame3$Protein_Length_aa, "aa\n")

# Compare our TP53 sequence with NCBI CDS

# Sequence length comparison
our_seq_length <- length(computational_orf_seq)
ncbi_seq_length <- length(ncbi_cds_seq)

cat("Our TP53 sequence length:", our_seq_length, "nt\n")
cat("NCBI CDS sequence length:", ncbi_seq_length, "nt\n")

# Compare nucleotide composition
our_nt_chars <- strsplit(as.character(computational_orf_seq), "")[[1]]
our_nt_counts <- table(our_nt_chars)
our_nt_percent <- (our_nt_counts / sum(our_nt_counts)) * 100

cat("\nOur TP53 nucleotide composition:\n")
print(our_nt_counts)
print(our_nt_percent)

cat("\nNCBI nucleotide composition:\n")
print(ncbi_nt_counts)
print(ncbi_nt_percent)

# Compare GC content
our_gc_content <- (our_nt_counts["G"] + our_nt_counts["C"]) /
  sum(our_nt_counts) * 100

cat("\nOur TP53 GC content:", our_gc_content, "%\n")
cat("NCBI GC content:", gc_content, "%\n")

# Check whether sequences are exactly identical
sequence_match <- identical(
  as.character(computational_orf_seq),
  as.character(ncbi_cds_seq)
)

cat("\nExact sequence match:", sequence_match, "\n")

# Compare computational ORF with NCBI-annotated CDS

# ORF coordinates
our_orf_start <- (longest_orf$start_codon - 1) * 3 + 2
our_orf_end <- (longest_orf$stop_codon - 1) * 3 + 2

# NCBI CDS coordinates
ncbi_cds_start <- 1
ncbi_cds_end <- ncbi_cds_length

cat("Our computational ORF coordinates:",
    our_orf_start, "-", our_orf_end, "\n")

cat("NCBI CDS coordinates:",
    ncbi_cds_start, "-", ncbi_cds_end, "\n")


# Compare ORF nucleotide sequence with NCBI CDS
our_orf_sequence <- substr(
  as.character(computational_orf_seq),
  our_orf_start,
  our_orf_end
)

ncbi_sequence <- as.character(ncbi_cds_seq)

sequence_identical <- identical(
  our_orf_sequence,
  ncbi_sequence
)

cat("\nORF nucleotide sequence identical to NCBI CDS:",
    sequence_identical, "\n")


# Translate our computational ORF
our_orf_protein <- Biostrings::translate(
  Biostrings::DNAString(our_orf_sequence)
)

# Translate NCBI CDS
ncbi_protein <- Biostrings::translate(
  Biostrings::DNAString(ncbi_sequence)
)

# Compare translated proteins
protein_identical <- identical(
  as.character(our_orf_protein),
  as.character(ncbi_protein)
)

cat("Translated protein identical:",
    protein_identical, "\n")

# Protein lengths excluding the terminal stop codon

our_orf_protein_length <- nchar(as.character(our_orf_protein))
ncbi_protein_length <- nchar(as.character(ncbi_protein))

if (endsWith(as.character(our_orf_protein), "*")) {
  our_orf_protein_length <- our_orf_protein_length - 1
}

if (endsWith(as.character(ncbi_protein), "*")) {
  ncbi_protein_length <- ncbi_protein_length - 1
}

cat("Our ORF protein length:",
    our_orf_protein_length, "aa\n")

cat("NCBI protein length:",
    ncbi_protein_length, "aa\n")

# Nucleotide composition plot

nucleotide_counts <- c(
  A = sum(orf_chars == "A"),
  T = sum(orf_chars == "T"),
  G = sum(orf_chars == "G"),
  C = sum(orf_chars == "C")
)

barplot(
  nucleotide_counts,
  main = "TP53 Nucleotide Composition",
  xlab = "Nucleotide",
  ylab = "Count"
)

nucleotide_percent <- nucleotide_counts / sum(nucleotide_counts) * 100

barplot(
  nucleotide_percent,
  main = "TP53 Nucleotide Composition (%)",
  xlab = "Nucleotide",
  ylab = "Percentage (%)"
)

# GC vs AT content visualization
gc_at_values <- c(
  GC = gc_content,
  AT = at_content
)
barplot(
  gc_at_values,
  main = "TP53 GC and AT Content",
  xlab = "Nucleotide Content",
  ylab = "Percentage (%)",
  ylim = c(0, 100)
)

# Professional nucleotide composition plot
library(ggplot2)

nucleotide_df <- data.frame(
  Nucleotide = names(nucleotide_counts),
  Count = as.numeric(nucleotide_counts)
)
ggplot(nucleotide_df, aes(x = Nucleotide, y = Count)) +
  geom_col() +
  labs(
    title = "TP53 Nucleotide Composition",
    x = "Nucleotide",
    y = "Count"
  ) +
  theme_minimal()

# Professional ORF protein length plot
orf_df <- data.frame(
  Reading_Frame = names(orf_protein_lengths),
  Protein_Length = as.numeric(orf_protein_lengths)
)

ggplot(orf_df, aes(x = Reading_Frame, y = Protein_Length)) +
  geom_col() +
  labs(
    title = "Longest ORF Protein Length by Reading Frame",
    x = "Reading Frame",
    y = "Protein Length (aa)"
  ) +
  theme_minimal()

# Compare ORF nucleotide and protein lengths
orf_comparison <- data.frame(
  Reading_Frame = c("Frame1", "Frame2", "Frame3"),
  Nucleotide_Length = c(1182, 204, 18),
  Protein_Length = c(393, 67, 5)
)

# Convert to long format
orf_long <- data.frame(
  Reading_Frame = rep(orf_comparison$Reading_Frame, 2),
  Length = c(
    orf_comparison$Nucleotide_Length,
    orf_comparison$Protein_Length
  ),
  Type = rep(c("Nucleotide (nt)", "Protein (aa)"), each = 3)
)

ggplot(orf_long, aes(x = Reading_Frame, y = Length, fill = Type)) +
  geom_col(position = "dodge") +
  labs(
    title = "ORF Nucleotide and Protein Length Comparison",
    x = "Reading Frame",
    y = "Length"
  ) +
  theme_minimal()

plot(
  gc_by_position,
  type = "l",
  main = "TP53 GC Content Across Sequence",
  xlab = "Sequence Position",
  ylab = "GC Content (%)"
)

# Nucleotide composition across codon positions
plot(
  1:3,
  gc_by_position,
  type = "o",
  pch = 16,
  ylim = c(0, 100),
  main = "TP53 Nucleotide Composition Across Codon Positions",
  xlab = "Codon Position",
  ylab = "GC Content (%)"
)

lines(
  1:3,
  100 - gc_by_position,
  type = "o",
  pch = 16
)

legend(
  "topright",
  legend = c("GC", "AT"),
  lty = 1,
  pch = 16
)

# Scatter plot: ORF nucleotide length vs protein length
library(ggplot2)

ggplot(orf_comparison, aes(x = Nucleotide_Length,
                           y = Protein_Length,
                           label = Reading_Frame)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  labs(
    title = "ORF Nucleotide Length vs Protein Length",
    x = "Nucleotide Length (nt)",
    y = "Protein Length (aa)"
  ) +
  theme_minimal()

# Heatmap: GC content across codon positions
library(ggplot2)

gc_heatmap <- data.frame(
  Codon_Position = c("Position 1", "Position 2", "Position 3"),
  GC_Content = as.numeric(gc_by_position)
)

ggplot(gc_heatmap, aes(x = Codon_Position, y = "TP53", fill = GC_Content)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = paste0(round(GC_Content, 1), "%")), size = 5) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue"
  ) +
  labs(
    title = "TP53 GC Content Across Codon Positions",
    x = "Codon Position",
    y = NULL,
    fill = "GC (%)"
  ) +
  theme_minimal()

# ============================================
# TP53 Nucleotide Composition: Our vs NCBI
# Heatmap
# ============================================

nucleotide_comparison <- data.frame(
  Nucleotide = c("A", "T", "G", "C"),
  Our_Data = as.numeric(our_nt_percent),
  NCBI = as.numeric(nt_percent)
)

# Convert to long format
library(tidyr)
library(ggplot2)

nucleotide_heatmap <- nucleotide_comparison %>%
  pivot_longer(
    cols = c(Our_Data, NCBI),
    names_to = "Dataset",
    values_to = "Percentage"
  )

# Heatmap
ggplot(nucleotide_heatmap,
       aes(x = Nucleotide,
           y = Dataset,
           fill = Percentage)) +
  
  geom_tile(color = "white", linewidth = 1) +
  
  geom_text(
    aes(label = paste0(round(Percentage, 1), "%")),
    size = 5
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "steelblue"
  ) +
  
  labs(
    title = "TP53 Nucleotide Composition: Our Data vs NCBI",
    x = "Nucleotide",
    y = NULL,
    fill = "Percentage (%)"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(size = 11)
  )

# INTERACTIVE VISUALIZATIONS USING PLOTLY
 # install.packages("plotly")
library(plotly)
library(ggplot2)

# 1. INTERACTIVE NUCLEOTIDE COMPOSITION
#    OUR DATA vs NCBI

nucleotide_comparison <- data.frame(
  Nucleotide = c("A", "C", "G", "T"),
  Our_Data = as.numeric(our_nt_percent),
  NCBI = as.numeric(nt_percent)
)

p1 <- ggplot(
  nucleotide_comparison,
  aes(x = Nucleotide)
) +
  geom_col(
    aes(
      y = Our_Data,
      text = paste(
        "Nucleotide:", Nucleotide,
        "<br>Our Data:", round(Our_Data, 1), "%"
      )
    ),
    position = position_dodge(width = 0.8),
    width = 0.35
  ) +
  geom_col(
    aes(
      y = NCBI,
      text = paste(
        "Nucleotide:", Nucleotide,
        "<br>NCBI:", round(NCBI, 1), "%"
      )
    ),
    position = position_nudge(x = 0.35),
    width = 0.35
  ) +
  labs(
    title = "TP53 Nucleotide Composition: Our Data vs NCBI",
    x = "Nucleotide",
    y = "Percentage (%)"
  ) +
  theme_minimal()

plotly_nucleotide <- ggplotly(
  p1,
  tooltip = "text"
)

plotly_nucleotide

# 2. INTERACTIVE GC CONTENT ACROSS CODON POSITIONS
gc_heatmap <- data.frame(
  Codon_Position = c(
    "Position 1",
    "Position 2",
    "Position 3"
  ),
  GC_Content = as.numeric(gc_by_position)
)

p2 <- ggplot(
  gc_heatmap,
  aes(
    x = Codon_Position,
    y = GC_Content,
    text = paste(
      "Codon Position:", Codon_Position,
      "<br>GC Content:", round(GC_Content, 1), "%"
    )
  )
) +
  geom_col(
    width = 0.6
  ) +
  labs(
    title = "TP53 GC Content Across Codon Positions",
    x = "Codon Position",
    y = "GC Content (%)"
  ) +
  theme_minimal()

plotly_gc <- ggplotly(
  p2,
  tooltip = "text"
)

plotly_gc

# 3. INTERACTIVE NUCLEOTIDE HEATMAP
#    OUR DATA vs NCBI

nucleotide_heatmap <- nucleotide_comparison

nucleotide_heatmap <- nucleotide_heatmap[
  rep(1:nrow(nucleotide_heatmap), 2),
]

nucleotide_heatmap$Dataset <- rep(
  c("Our Data", "NCBI"),
  each = 4
)

nucleotide_heatmap$Percentage <- c(
  nucleotide_comparison$Our_Data,
  nucleotide_comparison$NCBI
)

plotly_heatmap <- plot_ly(
  nucleotide_heatmap,
  x = ~Nucleotide,
  y = ~Dataset,
  z = ~Percentage,
  type = "heatmap",
  text = ~paste(
    "Dataset:", Dataset,
    "<br>Nucleotide:", Nucleotide,
    "<br>Percentage:", round(Percentage, 1), "%"
  ),
  hoverinfo = "text"
) %>%
  layout(
    title = "TP53 Nucleotide Composition: Our Data vs NCBI",
    xaxis = list(title = "Nucleotide"),
    yaxis = list(title = "")
  )

plotly_heatmap

# 4. INTERACTIVE ORF / PROTEIN LENGTH COMPARISON
if (exists("orf_df")) {
  
  p4 <- ggplot(
    orf_df,
    aes(
      x = Reading_Frame,
      y = Protein_Length,
      text = paste(
        "Reading Frame:", Reading_Frame,
        "<br>Protein Length:", Protein_Length,
        "aa"
      )
    )
  ) +
    geom_col() +
    labs(
      title = "TP53 Protein Length Across Reading Frames",
      x = "Reading Frame",
      y = "Protein Length (aa)"
    ) +
    theme_minimal()
  
  plotly_orf <- ggplotly(
    p4,
    tooltip = "text"
  )
  
  plotly_orf
}

# 5. INTERACTIVE AMINO ACID COMPOSITION
if (exists("tp53_aa_percent")) {
  
  # Get amino-acid names from the original count vector
  amino_acids <- names(tp53_aa_counts)
  
  # Make sure the number of labels matches the number of values
  if (is.null(amino_acids) || length(amino_acids) != length(tp53_aa_percent)) {
    
    amino_acids <- c(
      "A", "C", "D", "E", "F", "G", "H",
      "I", "K", "L", "M", "N", "P", "Q",
      "R", "S", "T", "V", "W", "Y", "*"
    )
  }
  
  aa_df <- data.frame(
    Amino_Acid = amino_acids,
    Percentage = as.numeric(tp53_aa_percent)
  )
  
  p5 <- ggplot(
    aa_df,
    aes(
      x = Amino_Acid,
      y = Percentage,
      text = paste(
        "Amino Acid:", Amino_Acid,
        "<br>Percentage:", round(Percentage, 2), "%"
      )
    )
  ) +
    geom_col(width = 0.7) +
    labs(
      title = "TP53 Amino Acid Composition",
      x = "Amino Acid",
      y = "Percentage (%)"
    ) +
    theme_minimal()
  
  plotly_aa <- ggplotly(
    p5,
    tooltip = "text"
  )
  
  plotly_aa
}

# 6. INTERACTIVE GC vs AT COMPOSITION
if (exists("gc_at")) {
  
  gc_at_df <- data.frame(
    Composition = c("GC", "AT"),
    Percentage = as.numeric(gc_at)
  )
  
  p6 <- ggplot(
    gc_at_df,
    aes(
      x = Composition,
      y = Percentage,
      text = paste(
        Composition,
        "Content:",
        round(Percentage, 1), "%"
      )
    )
  ) +
    geom_col(
      width = 0.6
    ) +
    labs(
      title = "TP53 GC vs AT Composition",
      x = "Base Composition",
      y = "Percentage (%)"
    ) +
    theme_minimal()
  
  plotly_gc_at <- ggplotly(
    p6,
    tooltip = "text"
  )
  
  plotly_gc_at
  
}

# END OF PLOTLY SECTION

 #install.packages("patchwork")
library(patchwork)

# PATCHWORK - COMBINE ALL ANALYSIS PLOTS

final_plot <- 
  (p1 | p2) /
  (p4 | p5) /
  (p6 | plot_spacer()) +
  plot_annotation(
    title = "TP53 DNA Sequence Analysis"
  ) &
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8)
  )

final_plot

# PATCHWORK - COMBINE ALL ANALYSIS PLOTS

library(patchwork)

final_plot <- (p1 | p2) /
  (p4 | p5) /
  p6 +
  plot_layout(
    heights = c(1, 1, 0.8)
  ) +
  plot_annotation(
    title = "TP53 DNA Sequence Analysis",
    theme = theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5
      )
    )
  )

final_plot

# Save final TP53 analysis results

sink("TP53_R_results.txt")

cat("TP53 DNA SEQUENCE ANALYSIS RESULTS\n")
cat("=================================\n\n")

cat("CDS length:", cds_length, "nt\n")
cat("Protein length:", protein_length, "aa\n")
cat("GC content:", gc_content, "%\n")
cat("AT content:", at_content, "%\n")
cat("Molecular weight:", our_mw, "Da\n")
cat("Theoretical pI:", our_pI, "\n\n")

cat("Longest ORFs:\n")
cat("Frame 1:", longest_orf_frame1$Length_nt, "nt /",
    longest_orf_frame1$Protein_Length_aa, "aa\n")
cat("Frame 2:", longest_orf_frame2$Length_nt, "nt /",
    longest_orf_frame2$Protein_Length_aa, "aa\n")
cat("Frame 3:", longest_orf_frame3$Length_nt, "nt /",
    longest_orf_frame3$Protein_Length_aa, "aa\n")

sink()

file.exists("TP53_R_results.txt")

dev.off()