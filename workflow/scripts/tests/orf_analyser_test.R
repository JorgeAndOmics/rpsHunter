################################################################################
# ORF-aware extraction and export from BLAST hits against a genome
# - Loads BLAST hits (parquet)
# - Infers strand and normalizes coordinates
# - Builds GRanges and harmonizes sequence IDs to FASTA headers
# - Extracts sequences honoring strand
# - Calls ORFs and keeps only slices containing ≥1 ORF
# - Exports: CSV (with ORF_COUNT and ORF_SUMMARY), GFF3, and FASTA
#
# Notes:
# - BLAST coordinates are 1-based inclusive.
# - Biostrings::getSeq respects strand: '-' yields reverse-complement.
# - ID harmonization handles prefixes like 'ref|...|', 'gb|...|' and drops .version.
################################################################################

# ==============================================================================
# 0) DEPENDENCIES
# ==============================================================================
options(warn = -1)  # silence warnings (flip to 0 while developing if desired)
suppressMessages({
  library(tidyverse)     # data wrangling, readr for CSV
  library(arrow)         # read_parquet
  library(GenomicRanges) # GRanges container and helpers
  library(plyranges)     # dplyr-like verbs for GRanges (not strictly required)
  library(rtracklayer)   # GFF3 export
  library(Biostrings)    # FASTA IO, getSeq, reverse-complement handling
  library(ORFik)         # ORF finding on DNAStringSet / GRanges
  library(IRanges)       # IRanges utils: elementNROWS, start/end on ranges
})

# ==============================================================================
# 1) CONFIGURATION
#    Edit these paths/parameters as needed or pass via commandArgs if preferred.
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)  # unused here but kept for compatibility
args.blast_table_output <- 'C:/Users/Lympha/Desktop/rpshunter_orf_test/tables/blast.parquet'
args.genome_fasta       <- 'V:/databases/local/blast_dbs_OLD/species/Hipposideros_larvatus.fa'
args.output_folder      <- 'C:/Users/Lympha/Desktop/rpshunter_orf_test/results'
args.min_orf_length     <- 200  # minimum ORF length in nucleotides

# Ensure output directory exists
dir.create(args.output_folder, showWarnings = FALSE, recursive = TRUE)

message("Starting BLAST genomic range and ORF analysis...")

# ==============================================================================
# 2) LOAD BLAST TABLE
#    Keep original column names as produced by BLAST to avoid confusion.
#    We:
#      - filter to the target Species
#      - drop rows without subject coordinates
#      - coerce S. Start / S. End to integer
# ==============================================================================
message("Loading BLAST data...")
d0 <- arrow::read_parquet(args.blast_table_output)

d1 <- d0 %>%
  filter(Species == "Hipposideros_larvatus") %>%                 # keep this species only
  filter(!is.na(`Subject ID`), !is.na(`S. Start`), !is.na(`S. End`)) %>%  # require coords
  mutate(
    `S. Start` = as.integer(`S. Start`),                         # BLAST subject start
    `S. End`   = as.integer(`S. End`)                            # BLAST subject end
  )

message(sprintf("Loaded %d hits after filtering.", nrow(d1)))

# ==============================================================================
# 3) STRAND INFERENCE + COORDINATE NORMALIZATION
#    BLAST uses 1-based inclusive coordinates. We:
#      - infer strand by ordering of S. Start and S. End
#      - define normalized [start, end] with start <= end for GRanges
# ==============================================================================
d2 <- d1 %>%
  mutate(
    strand = if_else(`S. Start` < `S. End`, "+", "-"),           # '+' if forward, '-' if reverse
    start  = pmin(`S. Start`, `S. End`),                         # normalized start (<= end)
    end    = pmax(`S. Start`, `S. End`)                          # normalized end
  )

# ==============================================================================
# 4) BUILD GRanges FROM THE DATA FRAME
#    - seqnames: 'Subject ID' must correspond to sequence IDs in the FASTA
#    - start/end: normalized genomic interval
#    - strand: '+' or '-'
#    - keep.extra.columns=TRUE keeps original BLAST columns in mcols(gr)
# ==============================================================================
gr <- makeGRangesFromDataFrame(
  d2,
  seqnames.field      = "Subject ID",
  start.field         = "start",
  end.field           = "end",
  strand.field        = "strand",
  keep.extra.columns  = TRUE
)

# ==============================================================================
# 5) READ GENOME FASTA AND HARMONIZE SEQUENCE IDENTIFIERS
#    FASTA headers and BLAST Subject IDs often differ:
#      - DB prefixes 'ref|ACC.VERSION|' or 'gb|...|' or 'emb|...|'
#      - 'chr' prefix in some builds
#      - version suffix '.N' at the end
#    We:
#      - normalize both GRanges seqlevels and FASTA names with the same function
#      - build a mapping from GRanges levels -> exact FASTA names
#      - rename GRanges levels where a match exists
#      - drop any remaining levels that are not present in the FASTA
# ==============================================================================
# Normalizer for IDs (vectorized)
normalize_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)                         # remove surrounding whitespace
  x <- sub("^[^|]*\\|", "", x)           # drop leading DB prefix like 'ref|'/'gb|'/'emb|'
  x <- sub("\\|$", "", x)                # drop trailing pipe if present
  x <- sub("^chr", "", x, ignore.case = TRUE) # drop 'chr' prefix (case-insensitive)
  x <- sub("\\.\\d+$", "", x)            # drop version suffix '.number' at end
  x
}

# Read the genome FASTA as DNAStringSet; names(genome) are the headers
genome <- readDNAStringSet(args.genome_fasta)

# Prepare normalized keys for mapping
gn_raw         <- names(genome)          # exact FASTA headers
gn_norm        <- normalize_ids(gn_raw)  # normalized FASTA headers
seq_lvls       <- seqlevels(gr)          # current GRanges levels (as is)
seq_lvls_norm  <- normalize_ids(seq_lvls)# normalized GRanges levels

# For each GRanges level (normalized), find its position in normalized FASTA names
# match() returns index in gn_norm or NA if not found
mapped <- gn_raw[ match(seq_lvls_norm, gn_norm) ]  # target exact FASTA names
ok     <- !is.na(mapped)                            # TRUE where a mapping exists

# Rename only the matched levels: setNames(values=new_names, names=old_levels)
if (any(ok)) {
  gr <- renameSeqlevels(gr, setNames(mapped[ok], seq_lvls[ok]))
}

# After renaming, keep only seqlevels that actually exist in the FASTA
keep <- intersect(seqlevels(gr), names(genome))  # safe set
drop <- setdiff(seqlevels(gr), keep)             # will be pruned
if (length(drop)) {
  message("Dropping ", length(drop), " unmapped seqlevels: ",
          paste(utils::head(drop, 12), collapse = ", "),
          if (length(drop) > 12) " ...")
}
gr <- keepSeqlevels(gr, keep, pruning.mode = "coarse")

# ==============================================================================
# 6) EXTRACT SEQUENCES FOR EACH RANGE (STRAND-AWARE)
#    Biostrings::getSeq(genome, gr):
#      - extracts [start, end] from the named sequence in 'genome'
#      - if strand == '-', returns reverse-complement automatically
# ==============================================================================
seqs <- getSeq(genome, gr)

# ==============================================================================
# 7) FIND ORFs IN EACH EXTRACTED SEQUENCE
#    ORFik::findORFs on a DNAStringSet returns an IRangesList:
#      - each element contains local ORF intervals for the corresponding slice
#      - with longestORF=TRUE, ORFik keeps the longest ORF per frame per slice
#      - minimumLength is in nucleotides
# ==============================================================================
orfs <- ORFik::findORFs(
  seqs,
  longestORF     = TRUE,
  minimumLength  = args.min_orf_length
)

# Some ORFik versions drop empty entries (i.e., no element for sequences with 0 ORFs).
# Build indices of sequences that DO have at least one ORF:
has_idx <- if (length(orfs) == length(seqs)) {
  # Full-length list was returned; count elements per slice
  which(IRanges::elementNROWS(orfs) > 0)
} else {
  # ORFik dropped empties; 'names(orfs)' carry the original 1-based indices as strings
  sort(unique(as.integer(names(orfs))))
}

# Filter GRanges and its backing data.frame rows to those with ≥1 ORF
gr_filt <- gr[has_idx]
d2_filt <- d2[has_idx, , drop = FALSE]

# Align the ORF list to the order of gr_filt for summarization
orfs_filt <- if (length(orfs) == length(seqs)) {
  orfs[has_idx]
} else {
  # Map chosen indices back into 'orfs' by name
  orfs[ match(has_idx, as.integer(names(orfs))) ]
}

# ==============================================================================
# 8) SUMMARIZE ORFs PER SLICE FOR CSV
#    We add:
#      - ORF_COUNT   : number of ORFs in that slice
#      - ORF_SUMMARY : "seqid:genomeStart-genomeEnd(strand):AAlen; ..."
#    Local ORF coordinates must be mapped back to genome coordinates:
#      '+' strand: genome = slice_start + local - 1
#      '-' strand: genome = slice_end   - local + 1  (mirrored)
# ==============================================================================
summarize_orfs <- function(ir, parent_gr) {
  # ir: IRanges of ORFs in local coordinates (relative to the extracted slice)
  # parent_gr: the GRanges range for that slice (genome coordinates + strand)
  if (length(ir) == 0) return(NA_character_)  # safeguard
  
  ps <- IRanges::start(ir)  # local starts
  pe <- IRanges::end(ir)    # local ends
  
  if (as.character(strand(parent_gr)) == "+") {
    # '+' strand: shift local to genome by adding (slice_start - 1)
    gs <- start(parent_gr) + ps - 1L
    ge <- start(parent_gr) + pe - 1L
  } else {
    # '-' strand: mirror around slice_end
    gs <- end(parent_gr) - pe + 1L
    ge <- end(parent_gr) - ps + 1L
  }
  
  # Rough amino-acid length: floor(nt_len / 3); ignores internal frameshifts
  aa <- floor((pe - ps + 1L) / 3L)
  
  # Build items like "seqid:123-456(-):111aa"
  items <- paste0(
    as.character(seqnames(parent_gr)), ":", gs, "-", ge,
    "(", as.character(strand(parent_gr)), "):", aa, "aa"
  )
  
  # Join multiple ORFs with ';'
  paste(items, collapse = ";")
}

# Add columns to the filtered data.frame that correspond to gr_filt rows
d2_filt$ORF_COUNT <- as.integer(IRanges::elementNROWS(orfs_filt))
d2_filt$ORF_SUMMARY <- vapply(
  seq_along(gr_filt),
  function(j) summarize_orfs(orfs_filt[[j]], gr_filt[j]),
  FUN.VALUE = character(1)
)

# ==============================================================================
# 9) EXPORT CSV OF FILTERED HITS + ORF COLUMNS
#    Contains: original BLAST columns + strand/start/end + ORF_COUNT/ORF_SUMMARY
# ==============================================================================
message("Writing ORF data table CSV...")
readr::write_csv(d2_filt, file.path(args.output_folder, "orf_data_table.csv"))

# ==============================================================================
# 10) EXPORT GFF3 OF FILTERED RANGES
#     Minimal region features with an ID that encodes locus + tag + strand.
#     Also pass through selected BLAST attributes if present.
#     To avoid huge files, drop very large text columns (e.g., 'Subject Sequence').
# ==============================================================================
message("Writing filtered ranges GFF3...")

# Drop heavy columns if present to keep GFF3 lean (optional but recommended)
if ("Subject Sequence" %in% colnames(mcols(gr_filt))) {
  mcols(gr_filt)$`Subject Sequence` <- NULL
}

# Set a simple feature type
mcols(gr_filt)$type <- "region"

# Build ID attribute: <seqid>_<start>-<end>|<Tag>|<strand>
tag_vec <- mcols(gr_filt)$Tag
if (is.null(tag_vec)) tag_vec <- NA_character_
tag_vec <- ifelse(is.na(tag_vec) | tag_vec == "", "NA", as.character(tag_vec))
mcols(gr_filt)$ID <- paste0(
  as.character(seqnames(gr_filt)), "_",
  start(gr_filt), "-", end(gr_filt), "|",
  tag_vec, "|",
  as.character(strand(gr_filt))
)

# Pass through a few useful BLAST columns if they exist in metadata
keep_attr <- intersect(
  c("Query ID","E-value","Bit Score","Pct Identity","Species","Tag"),
  colnames(mcols(gr))
)
if (length(keep_attr) > 0) {
  mcols(gr_filt)[keep_attr] <- mcols(gr)[has_idx, keep_attr, drop = FALSE]
}

# Write GFF3 file
rtracklayer::export(
  gr_filt,
  con    = file.path(args.output_folder, "filtered_ranges.gff3"),
  format = "GFF3"
)

# ==============================================================================
# 11) EXPORT FASTA OF FILTERED RANGES
#     getSeq(genome, gr_filt) respects strand; headers encode locus, tag, strand.
# ==============================================================================
message("Writing filtered regions FASTA...")

fa <- getSeq(genome, gr_filt)  # strand-aware extraction

# FASTA headers: <seqid>_<start>-<end>|<tag>|<strand>
hdr <- paste0(
  as.character(seqnames(gr_filt)), "_",
  start(gr_filt), "-", end(gr_filt), "|",
  tag_vec, "|",
  as.character(strand(gr_filt))
)
names(fa) <- hdr

writeXStringSet(fa, file.path(args.output_folder, "filtered_regions.fa"))

# ==============================================================================
# 12) DONE
# ==============================================================================
message("Done. Slices with ORFs: ", length(gr_filt), " of ", length(gr), " total.")
