############################################
# DADA2 MiFish (12S) pipeline - Wanda Bodnar
############################################

# MERGE OF JANE HALLAM'S AND NEOF's CODE 

## ========================================

## Configuration

# Define the source folder for raw FASTQ files and the destination for results
input_dir    <- "~/PhD/Jane data/Hallam_Thames_12S_raw/12S_raw"   # folder of FASTQs
output_dir   <- "~/PhD/Jane data/Hallam_Thames_12S_raw/12S_raw/Output"
blast_dir  <- "C:/blast_work" # Main directory for BLAST steps

# Specific path to the Cutadapt executable (external tool for primer removal)
cutadapt_exe <- "C:/Users/bodna/miniconda3/Scripts/cutadapt.exe" 

# Define the MiFish-U universal primers used for fish eDNA amplification
FWD <- "GTCGGTAAAACTCGTGCCAGC"   # MiFish-U forward primer
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # MiFish-U reverse primer

# Filtering/trimming: maxEE (expected errors) is the main quality gate
maxEE    <- c(5, 5)   # if a read is predicted to have more than 5 errors based on its quality score, it gets tossed out
truncLen <- c(160, 160) # defines every read to exactly 160 base pairs
truncQ   <- 2 # it cuts the read off the moment it encounters a base with a quality score of 2
minLen   <- 60
minOverlap <- 12 # when merging forward and reverse reads, they must overlap by at least 12 identical bases to be considered a valid pair    
maxMismatch <- 0      

# Reference databases (Miya) for assigning fish names to DNA sequences
train_tax_fasta     <- "C:/Users/bodna/Dropbox/Wanda Bodnar/Data/Practice/DADA2/references.12s.miya.dada.taxonomy.v258.fasta"
train_species_fasta <- "C:/Users/bodna/Dropbox/Wanda Bodnar/Data/Practice/DADA2/references.12s.miya.dada.species.v258.fasta"

# Optional extra metadata (join by Sample_number + Bio_rep)
metadata_excel <- "C:/Users/bodna/Dropbox/Wanda Bodnar/Data/Practice/DADA2/Thames_samples_metadata.xlsx"

# Plotting settings for visual summaries
rank_to_plot <- "Family"  # for taxonomic bar chart
topN_taxa    <- 15

set.seed(42)  # ensures the random nature of NMDS/DADA2 is reproducible, getting exact same results

## ============================

# Load required bioinformatics and data manipulation libraries
suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(ShortRead)
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(vegan)
  library(iNEXT)
  library(indicspecies)
  library(patchwork)
  library(stringr)
})

# Standardises file paths for Windows/Unix compatibility
normp <- function(p, mustWork = FALSE) normalizePath(path.expand(p), 
                                                     winslash = "/", mustWork = mustWork)

# Primer detection, generates all possible orientations of a primer to check for presence in reads
all_orients <- function(primer) {
  dna <- DNAString(primer)
  c(Forward   = as.character(dna),
    Complement= as.character(complement(dna)),
    Reverse   = as.character(reverse(dna)),
    RevComp   = as.character(reverseComplement(dna)))
}

# Counts occurrences of a primer in a specific FASTQ file
primer_hits <- function(primer, fn) {
  p <- DNAString(as.character(primer))
  s <- sread(readFastq(fn))
  sum(vcountPattern(p, s, fixed = FALSE) > 0)
}

# Extracts site, transect, and replicate info from complex filenames
parse_sample_ids <- function(ids) {
  core   <- sub("_S\\d+$", "", sub("-12S$", "", ids, ignore.case = TRUE))
  tokens <- strsplit(core, "-", fixed = TRUE)
  out <- lapply(tokens, function(v){
    n  <- length(v)
    br <- v[n]  # e.g. "1a", "2d"
    nums <- which(grepl("^[0-9]+$", v))
    sn <- if (length(nums)) as.integer(v[max(nums)]) else NA_integer_
    data.frame(Sample_number = sn,
               Bio_rep = br,
               Transect = sub("[a-z]$", "", br),
               Rep = sub("^[0-9]+", "", br),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# Formats the taxonomy table into the standard kingdom-to-species matrix
build_tax_mat <- function(seqs, taxa_df) {
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  mat <- matrix(NA_character_, nrow = length(seqs), ncol = length(ranks),
                dimnames = list(seqs, ranks))
  if (!is.null(taxa_df)) {
    keep <- intersect(rownames(taxa_df), seqs)
    cols <- intersect(colnames(taxa_df), ranks)
    if (length(keep) && length(cols)) mat[keep, cols] <- as.matrix(taxa_df[keep, cols, drop = FALSE])
  }
  mat
}

# Setup paths
input_dir  <- normp(input_dir,  mustWork = TRUE)
output_dir <- normp(output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Input:  ", input_dir,  "\n")
cat("Output: ", output_dir, "\n")

# Create matched lists of Forward (R1) and Reverse (R2) file paths
fnFs <- sort(list.files(input_dir, pattern = "_L001_R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "_L001_R2_001\\.fastq(\\.gz)?$", full.names = TRUE))

# Ensures every forward read has a matching reverse partner
length(fnFs)
length(fnRs) 
length(fnFs) == length(fnRs)

stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))
cat("Found", length(fnFs), "pairs\n")

# DADA2 requires reads with zero ambiguous 'N' bases for error learning
# Illumina sequencer uses fluorescent signals to determine if the base is A,C,G, or T
# N is as a placeholder for a base that couldn't be identified
filtN_dir  <- file.path(output_dir, "filtN"); dir.create(filtN_dir, showWarnings = FALSE, recursive = TRUE)

# It ensures that the forward and reverse reads are perfectly paired and that their names are "clean"
sample_id <- function(x) sub("_L00[12]_R[12]_001\\.fastq(\\.gz)?$", "", basename(x))
sidF <- sample_id(fnFs); sidR <- sample_id(fnRs); stopifnot(identical(sidF, sidR))

# Selection criteria for the Thames sites and 3 biological replicates
keep_sites      <- c(1, 2, 3, 4, 5, 6, 
                     15, 16, 17, 
                     21, 23, 24, 25, 26, 27, 
                     37, 38, 39) # set to NULL for all
keep_transects  <- c(1, 2, 3)          # <-- 1x and 3x; add 2 if you also want 2x
letters_regex   <- "[a-z]"          # replicates a–z (use "[abce]" to restrict)

# Constructing a pattern to subset only the specific sites/transects intended for analysis
pat <- paste0(
  "-(", paste(keep_sites, collapse="|"), ")-(",
  paste(keep_transects, collapse="|"), ")", letters_regex, "-12S_")

# Identifying the matches
idx <- grep(pat, sidF, ignore.case = TRUE)
if (length(idx) == 0) stop("No matches — check site numbers/transects/letters.")

fnFs <- fnFs[idx]
fnRs <- fnRs[idx]
sidF <- sidF[idx]

cat("Keeping", length(idx), "paired samples:\n")
print(sidF)

## Recreate matching output paths AFTER subsetting
filtN_dir  <- file.path(output_dir, "filtN"); dir.create(filtN_dir, TRUE, FALSE)
fnFs.filtN <- file.path(filtN_dir, basename(fnFs))
fnRs.filtN <- file.path(filtN_dir, basename(fnRs))

# Executes the N-base removal
out.filter <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Expect to keep >95% of the reads during the N-filtering step
stats <- as.data.frame(out.filter)
stats$pct_kept <- (stats$reads.out / stats$reads.in) * 100
print(stats)

# Diagnostic step to ensure primers are actually present before attempting to trim them
FWD.orients <- all_orients(FWD); REV.orients <- all_orients(REV)
cat("\nPrimer hits (first pair):\n")
print(rbind(
  FwdInR1 = sapply(FWD.orients, primer_hits, fn = fnFs.filtN[[1]]),
  FwdInR2 = sapply(FWD.orients, primer_hits, fn = fnRs.filtN[[1]]),
  RevInR1 = sapply(REV.orients, primer_hits, fn = fnFs.filtN[[1]]),
  RevInR2 = sapply(REV.orients, primer_hits, fn = fnRs.filtN[[1]])))

# Cutadapt removes the primer sequences; these are biological markers, not fish DNA
stopifnot(length(tryCatch(system2(cutadapt_exe, "--version", stdout = TRUE), error = function(e) "")) > 0)
cat("\nCutadapt:", paste(tryCatch(system2(cutadapt_exe, "--version", stdout = TRUE), error = function(e) ""), collapse=" "), "\n")

path.cut <- normp(file.path(output_dir, "cutadapt"))
dir.create(path.cut, showWarnings = FALSE, recursive = TRUE)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Helper function to run Cutadapt command-line tool from within R
run_cutadapt <- function(i, verbose = FALSE) {
  args <- c(
    "-g", FWD, "-a", REV.RC,      # R1 5' & 3'
    "-G", REV, "-A", FWD.RC,      # R2 5' & 3'
    "-n", "2",
    "--discard-untrimmed",        # Crucial: removes reads where primers weren't found
    "-m", as.character(minLen),   # Discards reads that become too short after trimming
    "-o", shQuote(fnFs.cut[i]),
    "-p", shQuote(fnRs.cut[i]),
    shQuote(fnFs.filtN[i]),
    shQuote(fnRs.filtN[i])
  )
  if (verbose) cat("\ncutadapt", paste(args, collapse=" "), "\n")
  system2(cutadapt_exe, args = args)
}

run_cutadapt(1, verbose = TRUE)
stopifnot(file.exists(fnFs.cut[1]), file.exists(fnRs.cut[1]))
for (i in seq_along(fnFs)) run_cutadapt(i)

# Re-filtering to apply the maxEE threshold and remove low-quality tail ends
cutFs <- sort(list.files(path.cut, pattern = "_L001_R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_L001_R2_001\\.fastq(\\.gz)?$", full.names = TRUE))
stopifnot(length(cutFs) == length(cutRs), length(cutFs) > 0)

# Inspect quality - visual check
# Green line is the mean quality score -> used to decide where to trim
# Orange line is the median quality score
# X axis is the position in the read (1 to 250 base pairs)
# Y axis is the quality score
plotQualityProfile(cutFs[1:1])
plotQualityProfile(cutRs[1:1])

filt_dir <- normp(file.path(output_dir, "filtered"))
dir.create(filt_dir, showWarnings = FALSE)

filtFs <- file.path(filt_dir, basename(cutFs))
filtRs <- file.path(filt_dir, basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs,
                     maxN = 0, maxEE = maxEE, truncQ = truncQ,
                     minLen = minLen, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# Convert the LATEST filtering results to a dataframe
stats_final <- as.data.frame(out)

# Calculate the percentage of reads that survived the Quality Gate
stats_final$pct_kept <- (stats_final$reads.out / stats_final$reads.in) * 100

# Print the results
# 80% or more, filtering settings are solid
print(stats_final)

# Core DADA2 process: learns the error profile of the specific sequencing run
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Check: The black line (model) should fit the black dots (observed data) perfectly.
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

samp_names <- sub("_L00[12]_R[12]_001\\.fastq(\\.gz)?$", "", basename(cutFs))
exists_idx <- file.exists(filtFs) & file.exists(filtRs)

# Denoising: Distinguishes biological sequences from sequencing errors
derepFs <- derepFastq(filtFs[exists_idx], verbose = TRUE)
names(derepFs) <- samp_names[exists_idx]

derepRs <- derepFastq(filtRs[exists_idx], verbose = TRUE)
names(derepRs) <- samp_names[exists_idx]

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge forward and reverse reads into a single full-length 12S fragment
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
seqtab  <- makeSequenceTable(mergers) # turns data into a spreadsheet

# Fragment length check
# Sequences around 170 bp are fish - BLAST to verify
# Sequences around 253 bp are human
# Low and very high sequences are chimeras
table(nchar(getSequences(seqtab)))

# Remove Chimeras: Artificial DNA molecules created during PCR
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

# Should be over 90%
cat("\nNon-chimeric proportion:", sum(seqtab.nochim)/sum(seqtab), "\n")

# Creates a table showing how many reads survived each step of the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- samp_names
track <- as.data.frame(track)
print(track)

# Assigns fish species names using the Naive Bayes classifier and Miya database
taxa <- NULL
if (!is.null(train_tax_fasta) && file.exists(path.expand(train_tax_fasta))) {
  taxa <- assignTaxonomy(seqtab.nochim, path.expand(train_tax_fasta), multithread = TRUE, verbose = TRUE)
  if (!is.null(train_species_fasta) && file.exists(path.expand(train_species_fasta))) {
    taxa <- addSpecies(taxa, path.expand(train_species_fasta))
  }
}

# Saves final abundance and taxonomy tables for ecological analysis
tables_dir <- normp(file.path(output_dir, "tables"))
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Stable ASV IDs (replaces long DNA strings with "ASV_1", "ASV_2", etc.)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0("ASV_", seq_along(asv_seqs))

# Export FASTA for BLAST
asv_dna <- DNAStringSet(asv_seqs)
names(asv_dna) <- asv_headers
writeXStringSet(asv_dna, file.path(blast_dir, "ASVs.fasta"))

# Export TSV count table - abundance table
tab_counts <- t(seqtab.nochim)
rownames(tab_counts) <- asv_headers

# Calculate total counts for each ASV across all samples
total_counts <- rowSums(tab_counts)

# Combine with sequence lengths
asv_summary <- data.frame(
  ASV = rownames(tab_counts),
  TotalAbundance = total_counts,
  Length = nchar(asv_seqs))

# Sort by abundance
asv_summary <- asv_summary[order(-asv_summary$TotalAbundance), ]

head(asv_summary, 10)

write.table(tab_counts, file.path(blast_dir, "DADA2_output_ASV_counts.tsv"), 
            sep="\t", quote=FALSE, col.names=NA)

cat("Success: ASVs.fasta and counts.tsv exported to C:/blast_work/\n")


############################################
# RUN BLAST (OUTSIDE R)
############################################
# In CMD or PowerShell:
# cd C:\blast_work
# "C:\Program Files\NCBI\blast\bin\makeblastdb.exe" -in ASVs.fasta -dbtype nucl -out ASVs_DB
# "C:\Program Files\NCBI\blast\bin\blastn.exe" -query ASVs.fasta -db ASVs_DB -out ASV_matches.txt -outfmt 6

# Load the BLAST results
blast_file <- "C:/blast_work/ASV_matches.txt"
blast_res <- read.table(blast_file, header=FALSE, stringsAsFactors=FALSE)
colnames(blast_res) <- c("query", "subject", "identity", "alignment_length", 
                         "mismatches", "gap_opens", "q_start", "q_end", 
                         "s_start", "s_end", "evalue", "bit_score")

# Filter for high-confidence hits (85% identity as per your protocol)
best_blast <- blast_res %>%
  filter(identity >= 85) %>%
  group_by(query) %>%
  slice_max(order_by = identity, n = 1) %>%
  select(query, identity, subject)

# Join with your DADA2 Taxonomy
tax_df <- as.data.frame(taxa)
tax_df$Sequence <- rownames(tax_df)
tax_df$ASV_ID <- paste0("ASV_", seq_len(nrow(tax_df)))

final_comparison <- tax_df %>%
  left_join(best_blast, by = c("ASV_ID" = "query"))

# Write this to a CSV so you can check it manually in Excel
write.csv(final_comparison, "C:/blast_work/Taxonomy_Verification.csv")

# Create a 'Verified' Species column
# If BLAST found a match, use that; otherwise, keep the DADA2 assignment
final_comparison <- final_comparison %>%
  mutate(Final_Species = ifelse(!is.na(subject), subject, Species))

# Prepare the matrix for Phyloseq
# Keep the standard ranks but replace species with verified ones
tax_mat <- final_comparison %>%
  select(Kingdom, Phylum, Class, Order, Family, Genus, Final_Species) %>%
  as.matrix()

# CRITICAL: Row names must be the DNA sequences to match the OTU table
rownames(tax_mat) <- final_comparison$Sequence
colnames(tax_mat)[7] <- "Species" # Rename back to standard 'Species'

# Extract the ASV sequences from seqtab.nochim
asv_seqs <- colnames(seqtab.nochim)

# Subset tax_mat to only include ASVs in seqtab.nochim
tax_mat <- tax_mat[rownames(tax_mat) %in% asv_seqs, , drop = FALSE]

# Ensure the order matches seqtab.nochim columns
tax_mat <- tax_mat[match(asv_seqs, rownames(tax_mat)), , drop = FALSE]

# Now assign unique ASV IDs
asv_headers <- paste0("ASV_", seq_along(asv_seqs))
colnames(seqtab.nochim) <- asv_headers
rownames(tax_mat) <- asv_headers

seq_names <- rownames(seqtab.nochim)
head(seq_names)

meta_from_seq <- tibble(
  Full_name = seq_names,
  Bio_rep = str_extract(seq_names, "(?<=-)(\\d+[a-z])(?=-12S)"),
  Sample_number = as.integer(str_extract(seq_names, "(?<=-JH-)(\\d+)(?=-\\d+[a-z]-12S)")))

head(meta_from_seq, 10)

meta_from_seq <- meta_from_seq %>%
  rename(Sample_number_int = Sample_number)  

map_tbl <- tibble::tribble(
  ~Sample_number_int, ~Season,  ~Location,
  1, "Spring", "Richmond",
  2,  "Spring", "Kew Bridge",     
  3, "Spring", "Chiswick",
  4, "Spring", "Battersea",
  5, "Spring", "Greenwich, Cutty Sark",
  6, "Spring", "Kew Bridge",
  
  15, "Fall", "Teddington Lock US", 
  16, "Fall", "Gravesend, beluga",
  17, "Fall", "Putney bridge", 
  
  21, "Fall", "Billingsgate",
  23, "Fall", "Greys, Thurrock",
  24, "Fall", "Vauxhall", 
  25, "Fall", "Fulham",
  26, "Fall", "Barking, Lower Roding",
  27, "Fall", "North Woolwich",
  
  37, "Fall", "Chiswick",
  38, "Fall", "Kew Bridge", 
  39, "Fall", "Richmond")

meta_full <- meta_from_seq %>%
  left_join(map_tbl, by = "Sample_number_int") %>%
  mutate(
    Season   = factor(Season, levels = c("Spring", "Fall")),
    Location = factor(Location, levels = c("Teddington Lock US", "Richmond", "Kew Bridge", "Chiswick", 
                                           "Putney bridge", "Fulham", "Vauxhall", "Battersea", 
                                           "Billingsgate", "Greenwich, Cutty Sark", "North Woolwich", "Barking, Lower Roding", 
                                           "Gravesend, beluga", "Greys, Thurrock")))

meta_full %>% 
  filter(is.na(Location)) %>% 
  select(Sample_number_int) %>% 
  unique()

# Should be 0
sum(is.na(meta_full$Location))

meta_df <- meta_full %>%
  as.data.frame() %>%
  tibble::column_to_rownames(var = "Full_name") 

all(rownames(seqtab.nochim) == rownames(meta_df))  

dim(seqtab.nochim)
otu_table(seqtab.nochim, taxa_are_rows = FALSE)
rowSums(seqtab.nochim)[1:10]
setdiff(rownames(seqtab.nochim), rownames(meta_df))  # should be character(0)
setdiff(rownames(meta_df), rownames(seqtab.nochim))  # should be character(0)

seqtab_ps <- t(seqtab.nochim)  # now rows = ASVs, cols = samples

all(rownames(tax_mat) == colnames(seqtab.nochim))  # should be TRUE

ps_final <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),  # samples = rows
  sample_data(meta_df),
  tax_table(tax_mat))  # add taxonomy)

# Save the complete object - this is your "Gold Standard" data file
saveRDS(ps_final, file.path(output_dir, "Thames_Fish_Final_Verified.rds"))

otu_mat <- as.data.frame(otu_table(ps_final))
head(otu_mat[1:5, 1:5])  # check rows = samples, columns = ASVs

otu_with_meta <- otu_mat %>%
  tibble::rownames_to_column("Full_name") %>%  
  left_join(meta_full %>% select(Full_name, Location, Season, Bio_rep), 
            by = "Full_name")

head(otu_with_meta)

rownames(tax_mat) <- colnames(seqtab.nochim)

tax_df <- as.data.frame(tax_table(ps_final))
tax_df$ASV_ID <- rownames(tax_df)
head(tax_df)

tax_df <- tax_df %>%
  mutate(
    Species_display = case_when(
      !is.na(Species) & Species != ASV_ID ~ Species,  # keep valid species
      !is.na(Genus) & Genus != "" ~ Genus,           # fallback to genus
      !is.na(Family) & Family != "" ~ Family,       # fallback to family
      TRUE ~ "Unknown"                               # fallback if nothing
    )
  ) %>%
  select(ASV_ID, Species_display)

otu_long <- otu_mat %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Full_name") %>%
  tidyr::pivot_longer(
    cols = starts_with("ASV_"),
    names_to = "ASV_ID",
    values_to = "count")

otu_long <- otu_long %>%
  left_join(tax_df, by = "ASV_ID") %>%  # now ASV_ID exists in tax_df
  left_join(meta_full %>% select(Full_name, Location, Season), by = "Full_name")

# Final table
otu_species <- otu_long %>%
  group_by(Full_name, Location, Season, Species_display) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  pivot_wider(
    names_from = Species_display,
    values_from = count,
    values_fill = 0)

head(otu_species)


#############################################


# --- iNEXT ---

# Combining counts and metadata, then summing replicates by Location/Season
iNEXT_input_grouped <- otu_species %>%
  group_by(Location, Season) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  mutate(Assemblage_ID = paste(Location, Season, sep = " @ "))

# Transpose for iNEXT (species as rows, site-season as columns)
iNEXT_data <- t(iNEXT_input_grouped %>% select(-Location, -Season, -Assemblage_ID))
colnames(iNEXT_data) <- iNEXT_input_grouped$Assemblage_ID

# Run iNEXT Analysis (datatype = "abundance" is used because of the raw ASV counts)
iNEXT_out <- iNEXT(iNEXT_data, q = c(0, 1, 2), datatype = "abundance")

# Clean and filter data for plotting
inext_final_points <- fortify(iNEXT_out, type = 1) %>%
  separate(Assemblage, into = c("Location", "Season"), sep = " @ ") %>%
  filter(Method == "Observed") %>%
  mutate(
    Location = factor(Location, levels = location_levels),
    Order.q = factor(Order.q, levels = c("0", "1", "2"), labels = hill_labels),
    Season = factor(Season, levels = season_levels))

# Plot
ggplot(inext_final_points, aes(x = Location, y = y, color = Location)) +
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin = y.lwr, ymax = y.upr), width = 0.3, linewidth = 1) +
  facet_grid(rows = vars(Order.q), cols = vars(Season), scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.text = element_text(face = "bold")) +
  scale_x_discrete(drop = FALSE) + 
  labs(
    title = "Fish diversity across the Thames",
    y = "Hill diversity (q = 0, 1, 2)",
    x = "Location",
    caption = "Confidence intervals (95%) are narrow due to high sequencing coverage.")


# --- NMDS ---
metadata <- otu_species[, 1:3]
species_matrix <- otu_species[, 4:ncol(otu_species)]
species_matrix <- species_matrix[, !grepl("Unknown", colnames(species_matrix), ignore.case = TRUE)]
row_totals <- rowSums(species_matrix)

species_final <- species_matrix[row_totals > 0, ]
metadata_final <- otu_species[row_totals > 0, 1:3]

# Bray
spec_hel <- decostand(species_final, method = "hellinger")
nmds_bray <- metaMDS(spec_hel, distance = "bray", k = 2, trymax = 100)

scores_bray <- as.data.frame(scores(nmds_bray, "sites"))
scores_bray$Location <- metadata_final$Location
scores_bray$Season <- metadata_final$Season

plot_bray <- ggplot(scores_bray, aes(x = NMDS1, y = NMDS2, color = Location, shape = Season)) +
  geom_point(size = 3) + 
  stat_ellipse() + 
  theme_bw() +
  labs(title = "Bray-Curtis (Hellinger)", subtitle = paste("Stress:", round(nmds_bray$stress, 3)))

adonis2(spec_hel ~ Location * Season, data = metadata_final, method = "bray")

plot_bray

# Jaccard
nmds_jaccard <- metaMDS(species_final, distance = "jaccard", binary = TRUE, k = 2, trymax = 100)

scores_jac <- as.data.frame(scores(nmds_jaccard, "sites"))
scores_jac$Location <- metadata_final$Location
scores_jac$Season <- metadata_final$Season

plot_jac <- ggplot(scores_jac, aes(x = NMDS1, y = NMDS2, color = Location, shape = Season)) +
  geom_point(size = 3) + 
  stat_ellipse() + 
  theme_bw() +
  labs(title = "Jaccard (Presence-Absence)", subtitle = paste("Stress:", round(nmds_jaccard$stress, 3)))

adonis2(species_final ~ Location * Season, data = metadata_final, method = "jaccard", binary = TRUE)

plot_jac

# PCoA
dist_matrix <- vegdist(spec_hel, method = "bray")

# Run the PCoA
pcoa_res <- wcmdscale(dist_matrix, eig = TRUE)

# Calculate the percentage explained by the first two axes
pc1_exp <- round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1)
pc2_exp <- round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1)

# Extract scores for plotting
scores_pcoa <- as.data.frame(pcoa_res$points)
colnames(scores_pcoa) <- c("PCoA1", "PCoA2")
scores_pcoa$Location <- metadata_final$Location
scores_pcoa$Season <- metadata_final$Season

# Plot
plot_pcoa <- ggplot(scores_pcoa, aes(x = PCoA1, y = PCoA2, color = Location, shape = Season)) +
  geom_point(size = 3) +
  stat_ellipse() +
  theme_bw() +
  labs(title = "PCoA (Bray-Curtis)", 
       x = paste0("PCoA1 (", pc1_exp, "%)"), 
       y = paste0("PCoA2 (", pc2_exp, "%)"))

plot_pcoa


# --- Ind ---
inv <- multipatt(spec_hel, metadata_final$Location, func = "IndVal.g", control = how(nperm = 999))
summary(inv)

heatmap_data <- cbind(metadata_final, spec_hel)

heatmap_avg <- heatmap_data %>%
  group_by(Location, Season) %>%
  summarise(across(where(is.numeric), mean)) %>%
  pivot_longer(cols = -c(Location, Season), names_to = "Species", values_to = "Abundance") %>%
  mutate(Species = reorder(Species, Abundance, sum))

# Bray
ggplot(heatmap_avg, aes(x = Location, y = Species, fill = Abundance)) +
  geom_tile(color = "white") +
  facet_wrap(~Season) +  
  scale_fill_gradient(low = "white", high = "darkred") +
  theme_minimal() +
  labs(title = "Indicator species heatmap",
       subtitle = "Mean hellinger-transformed eDNA abundance",
       fill = "Relative\nabundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


species_pa <- ifelse(species_final > 0, 1, 0)
heatmap_data_jac <- cbind(metadata_final, species_pa)

heatmap_avg_jac <- heatmap_data_jac %>%
  group_by(Location, Season) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  pivot_longer(cols = -c(Location, Season), names_to = "Species", values_to = "Frequency") %>%
  mutate(Species = reorder(Species, Frequency, sum))

# Jaccard
ggplot(heatmap_avg_jac, aes(x = Location, y = Species, fill = Frequency)) +
  geom_tile(color = "white") +
  facet_wrap(~Season) +  
  scale_fill_gradient(low = "white", high = "darkgreen") + 
  theme_minimal() +
  labs(title = "Indicator species heatmap",
       subtitle = "Proportion of samples where species was detected (0.0 to 1.0) (Jaccard-style)",
       fill = "Detection\nfrequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

