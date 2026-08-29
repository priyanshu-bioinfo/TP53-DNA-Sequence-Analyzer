# Human TP53 DNA Sequence Analyzer
# Reproducible R/RStudio component
#
# Run from the repository root in RStudio or with:
#   Rscript Code/R/DNA_analysis.R
#
# Required packages:
#   Biostrings, Peptides, ggplot2, patchwork

suppressPackageStartupMessages({
  library(Biostrings)
  library(Peptides)
  library(ggplot2)
  library(patchwork)
})

# -------------------------------------------------------------------------
# 1. Locate repository root and input/output files
# -------------------------------------------------------------------------

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "Data", "TP53.fa")) &&
        file.exists(file.path(current, "Data", "TP53_CDS.fa"))) {
      return(current)
    }

    parent <- dirname(current)
    if (parent == current) {
      stop("Could not find project root. Run this script from the repository.")
    }
    current <- parent
  }
}

ROOT <- find_project_root()
DATA <- file.path(ROOT, "Data")
RESULTS <- file.path(ROOT, "Results", "R")
PLOTS <- file.path(ROOT, "Plots", "R")

dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 2. Read the project sequences
# -------------------------------------------------------------------------

transcript <- readDNAStringSet(file.path(DATA, "TP53.fa"))
cds <- readDNAStringSet(file.path(DATA, "TP53_CDS.fa"))
protein_file <- readAAStringSet(file.path(DATA, "TP53_protein.fa"))
uniprot_reference <- readAAStringSet(file.path(DATA, "TP53_uniprot_reference.fa"))

if (length(transcript) != 1 || length(cds) != 1 ||
    length(protein_file) != 1 || length(uniprot_reference) != 1) {
  stop("Each reference FASTA file must contain exactly one sequence.")
}

transcript_seq <- toupper(as.character(transcript[[1]]))
cds_seq <- toupper(as.character(cds[[1]]))
protein_with_stop <- toupper(as.character(protein_file[[1]]))
protein_seq <- sub("\\*$", "", protein_with_stop)
reference_protein <- toupper(as.character(uniprot_reference[[1]]))
reference_protein <- sub("\\*$", "", reference_protein)

# Validate DNA alphabet.
if (grepl("[^ACGT]", transcript_seq) || grepl("[^ACGT]", cds_seq)) {
  stop("DNA input contains characters outside A/C/G/T.")
}

# -------------------------------------------------------------------------
# 3. Basic nucleotide analysis
# -------------------------------------------------------------------------

transcript_chars <- strsplit(transcript_seq, "", fixed = TRUE)[[1]]
cds_chars <- strsplit(cds_seq, "", fixed = TRUE)[[1]]

transcript_counts <- table(factor(transcript_chars, levels = c("A", "C", "G", "T")))
cds_counts <- table(factor(cds_chars, levels = c("A", "C", "G", "T")))

transcript_percent <- transcript_counts / sum(transcript_counts) * 100
cds_percent <- cds_counts / sum(cds_counts) * 100

transcript_gc <- (transcript_counts["G"] + transcript_counts["C"]) /
  sum(transcript_counts) * 100
transcript_at <- (transcript_counts["A"] + transcript_counts["T"]) /
  sum(transcript_counts) * 100

cds_gc <- (cds_counts["G"] + cds_counts["C"]) /
  sum(cds_counts) * 100
cds_at <- (cds_counts["A"] + cds_counts["T"]) /
  sum(cds_counts) * 100

# -------------------------------------------------------------------------
# 4. Reference validation
# -------------------------------------------------------------------------

# NCBI NM_000546.6 CDS coordinates recorded for this project are 143-1324.
annotated_cds_from_transcript <- substr(transcript_seq, 143, 1324)

cds_matches_transcript <- identical(cds_seq, annotated_cds_from_transcript)

translated_cds <- as.character(translate(DNAString(cds_seq)))
translated_cds_no_stop <- sub("\\*$", "", translated_cds)

protein_file_matches <- identical(protein_seq, reference_protein)
translated_cds_matches <- identical(translated_cds_no_stop, reference_protein)

# -------------------------------------------------------------------------
# 5. ORF detection in the full transcript
# -------------------------------------------------------------------------

find_orfs <- function(sequence) {
  sequence <- toupper(sequence)
  results <- list()

  for (frame_offset in 0:2) {
    frame_seq <- substr(
      sequence,
      frame_offset + 1,
      nchar(sequence)
    )

    usable_length <- nchar(frame_seq) - (nchar(frame_seq) %% 3)
    if (usable_length < 3) next

    frame_seq <- substr(frame_seq, 1, usable_length)
    codons <- substring(
      frame_seq,
      seq(1, usable_length - 2, by = 3),
      seq(3, usable_length, by = 3)
    )

    start_positions <- which(codons == "ATG")

    for (start_idx in start_positions) {
      stop_positions <- which(
        codons %in% c("TAA", "TAG", "TGA") &
          seq_along(codons) > start_idx
      )

      if (length(stop_positions) == 0) next

      stop_idx <- stop_positions[1]

      start_nt <- frame_offset + (start_idx - 1) * 3 + 1
      stop_nt <- frame_offset + stop_idx * 3
      length_nt <- stop_nt - start_nt + 1
      protein_length <- stop_idx - start_idx

      results[[length(results) + 1]] <- data.frame(
        Frame = paste0("+", frame_offset + 1),
        Start_Codon_Index = start_idx,
        Stop_Codon_Index = stop_idx,
        Start_nt = start_nt,
        Stop_nt = stop_nt,
        Stop_Codon = codons[stop_idx],
        Length_nt = length_nt,
        Protein_Length_aa = protein_length,
        stringsAsFactors = FALSE
      )

      # First in-frame stop after each ATG defines the candidate ORF.
    }
  }

  if (length(results) == 0) {
    return(data.frame(
      Frame = character(),
      Start_Codon_Index = integer(),
      Stop_Codon_Index = integer(),
      Start_nt = integer(),
      Stop_nt = integer(),
      Stop_Codon = character(),
      Length_nt = integer(),
      Protein_Length_aa = integer(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, results)
}

transcript_orfs <- find_orfs(transcript_seq)

if (nrow(transcript_orfs) == 0) {
  stop("No complete ORFs found in the TP53 transcript.")
}

longest_orf_by_frame <- do.call(
  rbind,
  lapply(split(transcript_orfs, transcript_orfs$Frame), function(x) {
    x[which.max(x$Length_nt), , drop = FALSE]
  })
)

longest_orf_overall <- transcript_orfs[
  which.max(transcript_orfs$Length_nt),
  , drop = FALSE
]

# The annotated TP53 CDS should correspond to transcript frame +2.
expected_orf_match <- with(
  longest_orf_overall,
  Start_nt == 143 &&
    Stop_nt == 1324 &&
    Length_nt == 1182
)

# -------------------------------------------------------------------------
# 6. Reading-frame analysis within the annotated CDS
# -------------------------------------------------------------------------

cds_frame_results <- lapply(0:2, function(offset) {
  frame_seq <- substr(cds_seq, offset + 1, nchar(cds_seq))
  usable_length <- nchar(frame_seq) - (nchar(frame_seq) %% 3)
  frame_seq <- substr(frame_seq, 1, usable_length)

  codons <- substring(
    frame_seq,
    seq(1, usable_length - 2, by = 3),
    seq(3, usable_length, by = 3)
  )

  start_positions <- which(codons == "ATG")
  rows <- list()

  for (start_idx in start_positions) {
    stops <- which(
      codons %in% c("TAA", "TAG", "TGA") &
        seq_along(codons) > start_idx
    )

    if (length(stops) == 0) next

    stop_idx <- stops[1]
    rows[[length(rows) + 1]] <- data.frame(
      Reading_Frame = paste0("Frame ", offset + 1),
      Start_Codon_Index = start_idx,
      Stop_Codon_Index = stop_idx,
      Length_nt = (stop_idx - start_idx + 1) * 3,
      Protein_Length_aa = stop_idx - start_idx,
      Stop_Codon = codons[stop_idx],
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) {
    return(data.frame(
      Reading_Frame = paste0("Frame ", offset + 1),
      Start_Codon_Index = integer(),
      Stop_Codon_Index = integer(),
      Length_nt = integer(),
      Protein_Length_aa = integer(),
      Stop_Codon = character(),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
})

cds_orfs <- do.call(rbind, cds_frame_results)

longest_cds_orfs <- do.call(
  rbind,
  lapply(split(cds_orfs, cds_orfs$Reading_Frame), function(x) {
    x[which.max(x$Length_nt), , drop = FALSE]
  })
)

# -------------------------------------------------------------------------
# 7. Codon and amino-acid analysis
# -------------------------------------------------------------------------

if (nchar(cds_seq) %% 3 != 0) {
  stop("CDS length is not divisible by 3.")
}

codons <- substring(
  cds_seq,
  seq(1, nchar(cds_seq) - 2, by = 3),
  seq(3, nchar(cds_seq), by = 3)
)

codon_counts <- table(codons)
sense_codons <- codons[!codons %in% c("TAA", "TAG", "TGA")]
sense_codon_counts <- table(sense_codons)
codon_percent <- prop.table(sense_codon_counts) * 100

aa_chars <- strsplit(reference_protein, "", fixed = TRUE)[[1]]
aa_counts <- table(factor(
  aa_chars,
  levels = strsplit("ACDEFGHIKLMNPQRSTVWY", "", fixed = TRUE)[[1]]
))
aa_percent <- aa_counts / sum(aa_counts) * 100

# GC content by codon position (sense codons only).
sense_matrix <- do.call(rbind, strsplit(sense_codons, "", fixed = TRUE))
gc_by_position <- c(
  "Position 1" = mean(sense_matrix[, 1] %in% c("G", "C")) * 100,
  "Position 2" = mean(sense_matrix[, 2] %in% c("G", "C")) * 100,
  "Position 3" = mean(sense_matrix[, 3] %in% c("G", "C")) * 100
)

# -------------------------------------------------------------------------
# 8. Protein physicochemical properties
# -------------------------------------------------------------------------

protein_length <- nchar(reference_protein)
protein_mw <- mw(reference_protein)
protein_pI <- pI(reference_protein)

# -------------------------------------------------------------------------
# 9. Plots
# -------------------------------------------------------------------------

p1_data <- data.frame(
  Nucleotide = names(cds_percent),
  Percentage = as.numeric(cds_percent)
)

p1 <- ggplot(p1_data, aes(Nucleotide, Percentage)) +
  geom_col() +
  labs(
    title = "TP53 CDS Nucleotide Composition",
    x = "Nucleotide",
    y = "Percentage (%)"
  ) +
  theme_minimal()

p2_data <- data.frame(
  Codon_Position = names(gc_by_position),
  GC_Content = as.numeric(gc_by_position)
)

p2 <- ggplot(p2_data, aes(Codon_Position, GC_Content)) +
  geom_col() +
  labs(
    title = "TP53 CDS GC Content Across Codon Positions",
    x = "Codon Position",
    y = "GC Content (%)"
  ) +
  theme_minimal()

p4_data <- longest_cds_orfs[, c(
  "Reading_Frame", "Protein_Length_aa"
)]
p4 <- ggplot(
  p4_data,
  aes(Reading_Frame, Protein_Length_aa)
) +
  geom_col() +
  labs(
    title = "TP53 Protein Length Across CDS Reading Frames",
    x = "Reading Frame",
    y = "Protein Length (aa)"
  ) +
  theme_minimal()

p5_data <- data.frame(
  Amino_Acid = names(aa_percent),
  Percentage = as.numeric(aa_percent)
)
p5 <- ggplot(p5_data, aes(Amino_Acid, Percentage)) +
  geom_col() +
  labs(
    title = "TP53 Amino Acid Composition",
    x = "Amino Acid",
    y = "Percentage (%)"
  ) +
  theme_minimal()

p6_data <- data.frame(
  Composition = c("GC", "AT"),
  Percentage = c(cds_gc, cds_at)
)
p6 <- ggplot(p6_data, aes(Composition, Percentage)) +
  geom_col() +
  labs(
    title = "TP53 CDS GC vs AT Composition",
    x = "Base Composition",
    y = "Percentage (%)"
  ) +
  ylim(0, 100) +
  theme_minimal()

# Save only after every plot object has been created.
ggsave(file.path(PLOTS, "p1.png"), p1, width = 10, height = 7, dpi = 300)
ggsave(file.path(PLOTS, "p2.png"), p2, width = 10, height = 7, dpi = 300)
ggsave(file.path(PLOTS, "p4.png"), p4, width = 10, height = 7, dpi = 300)
ggsave(file.path(PLOTS, "p5.png"), p5, width = 10, height = 7, dpi = 300)
ggsave(file.path(PLOTS, "p6.png"), p6, width = 10, height = 7, dpi = 300)

final_plot <- (p1 | p2) /
  (p4 | p5) /
  p6 +
  plot_annotation(title = "TP53 DNA Sequence Analysis")

ggsave(
  file.path(PLOTS, "final_plot.png"),
  final_plot,
  width = 12,
  height = 10,
  dpi = 300
)

# -------------------------------------------------------------------------
# 10. Save reproducible R results
# -------------------------------------------------------------------------

result_file <- file.path(RESULTS, "TP53_R_results.txt")

sink(result_file)
cat("TP53 DNA SEQUENCE ANALYSIS — R/RSTUDIO\n")
cat("======================================\n\n")

cat("INPUT / REFERENCE VALIDATION\n")
cat("----------------------------\n")
cat("Transcript ID:", names(transcript), "\n")
cat("Transcript length:", nchar(transcript_seq), "nt\n")
cat("CDS length:", nchar(cds_seq), "nt\n")
cat("CDS matches transcript positions 143-1324:", cds_matches_transcript, "\n")
cat("Translated CDS length:", nchar(translated_cds_no_stop), "aa\n")
cat("Translated CDS matches UniProt P04637:", translated_cds_matches, "\n")
cat("Protein file matches UniProt P04637:", protein_file_matches, "\n\n")

cat("NUCLEOTIDE COMPOSITION\n")
cat("----------------------\n")
for (base in c("A", "C", "G", "T")) {
  cat(
    base, ":", as.integer(cds_counts[base]),
    sprintf(" (%.2f%%)", cds_percent[base]), "\n"
  )
}
cat(sprintf("CDS GC content: %.5f%%\n", cds_gc))
cat(sprintf("CDS AT content: %.5f%%\n\n", cds_at))

cat("FULL-TRANSCRIPT ORF ANALYSIS\n")
cat("----------------------------\n")
for (frame_name in c("+1", "+2", "+3")) {
  x <- longest_orf_by_frame[longest_orf_by_frame$Frame == frame_name, ]
  if (nrow(x) == 1) {
    cat(
      "Frame", frame_name, ":",
      x$Start_nt, "-", x$Stop_nt, "nt;",
      x$Length_nt, "nt /",
      x$Protein_Length_aa, "aa;",
      "stop =", x$Stop_Codon, "\n"
    )
  }
}
cat("Longest ORF matches annotated TP53 CDS coordinates:", expected_orf_match, "\n\n")

cat("ANNOTATED CDS READING-FRAME ANALYSIS\n")
cat("------------------------------------\n")
for (i in seq_len(nrow(longest_cds_orfs))) {
  x <- longest_cds_orfs[i, ]
  cat(
    x$Reading_Frame, ":",
    x$Length_nt, "nt /",
    x$Protein_Length_aa, "aa;",
    "stop =", x$Stop_Codon, "\n"
  )
}
cat("\n")

cat("PROTEIN PROPERTIES — PEPTIDES PACKAGE\n")
cat("--------------------------------------\n")
cat("Protein length:", protein_length, "aa\n")
cat(sprintf("Molecular weight: %.2f Da\n", protein_mw))
cat(sprintf("Theoretical pI: %.6f\n\n", protein_pI))

cat("GC CONTENT BY CODON POSITION\n")
cat("----------------------------\n")
for (i in seq_along(gc_by_position)) {
  cat(
    names(gc_by_position)[i], ":",
    sprintf("%.2f%%", gc_by_position[i]), "\n"
  )
}
cat("\n")

cat("TOP 10 SENSE CODONS\n")
cat("-------------------\n")
top_codons <- sort(sense_codon_counts, decreasing = TRUE)
for (i in seq_len(min(10, length(top_codons)))) {
  cat(
    names(top_codons)[i], ":",
    as.integer(top_codons[i]), "\n"
  )
}
cat("\n")

cat("END OF R ANALYSIS\n")
sink()

cat("R analysis completed successfully.\n")
cat("Results:", result_file, "\n")
cat("Plots:", PLOTS, "\n")
