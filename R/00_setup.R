# packages, folders, helpers. run from the project root
stopifnot("open R in the project root: the folder that contains R/" = dir.exists("R"))

cran_pkgs <- c("tidyverse", "ggrepel")
bioc_pkgs <- c("GenomicRanges", "GenomicFeatures", "VariantAnnotation", "Biostrings",
               "TxDb.Hsapiens.UCSC.hg38.knownGene", "BSgenome.Hsapiens.UCSC.hg38")
have <- function(p) requireNamespace(p, quietly = TRUE)

miss <- cran_pkgs[!vapply(cran_pkgs, have, logical(1))]
if (length(miss)) install.packages(miss)
miss <- bioc_pkgs[!vapply(bioc_pkgs, have, logical(1))]
if (length(miss)) {
  if (!have("BiocManager")) install.packages("BiocManager")
  BiocManager::install(miss, ask = FALSE, update = FALSE)
}

# bioconductor first, tidyverse last so dplyr wins name clashes (select, rename, slice)
suppressPackageStartupMessages({
  library(GenomicRanges); library(GenomicFeatures); library(VariantAnnotation)
  library(Biostrings)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(BSgenome.Hsapiens.UCSC.hg38)
  library(tidyverse); library(ggrepel)
})

# folders and empty working files (written later, together)
for (d in c("data/raw", "docs", "notes", "results", "figures", "reports"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (f in c("docs/analysis_plan.md", "docs/data_dictionary.md",
            "docs/deviations_log.md", "notes/log.md"))
  if (!file.exists(f)) file.create(f)

# domain template: fill start/end from your uniprot file, delete or add rows as needed
if (!file.exists("data/domains_IL11RA.csv")) {
  tibble(feature = c("signal_peptide", "Ig-like_D1", "FNIII_D2", "FNIII_D3", "WSXWS_motif"),
         start = NA_integer_, end = NA_integer_,
         source = "UniProt Q14626", date_checked = NA_character_) |>
    write_csv("data/domains_IL11RA.csv", na = "")
}

# ---- helpers ----
locate <- function(pattern, dir = "data") {
  hits <- list.files(dir, pattern = pattern, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (!length(hits))
    stop("no file matching '", pattern, "' under ", dir, "/. files there now:\n  ",
         paste(list.files(dir, recursive = TRUE), collapse = "\n  "), call. = FALSE)
  if (length(hits) > 1) {
    hits <- hits[order(file.mtime(hits), decreasing = TRUE)]
    message("several matches, using the newest: ", hits[1])
  }
  hits[1]
}

check_cols <- function(df, cols, label) {
  miss <- setdiff(cols, names(df))
  if (length(miss))
    stop(label, ": missing column(s): ", paste(miss, collapse = ", "),
         "\nfound: ", paste(names(df), collapse = " | "), call. = FALSE)
}

AA3 <- c(A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", Q = "Gln", E = "Glu",
         G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys", M = "Met", F = "Phe",
         P = "Pro", S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val", "*" = "Ter")
AA1 <- setNames(names(AA3), AA3)

norm_p <- function(x) str_replace_all(str_remove_all(x, "[()\\s]"), "\\*", "Ter")

# alphamissense table (columns read by position, header names have odd characters)
read_am <- function() {
  f  <- locate("alphamissense.*\\.(zip|tsv)$")
  am <- read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE)
  stopifnot(ncol(am) >= 11, str_detect(am[[5]][1], "^p\\."))
  names(am)[1:11] <- c("uniprot", "entry", "gene", "enst", "protein_variant", "is_snv",
                       "aa_ref", "position", "aa_alt", "am_score", "am_class")
  am |> select(1:11) |> mutate(position = as.integer(position), am_score = as.numeric(am_score))
}

# uniprot Q14626 reference residues, index = position. uses a fasta in data/ if present,
# otherwise rebuilds the sequence from the alphamissense table (which is built on Q14626)
uniprot_ref <- function(am) {
  fa <- list.files("data", "\\.(fa|fasta)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(fa)) {
    s <- strsplit(as.character(readAAStringSet(fa[1])[[1]]), "")[[1]]
    message("reference: ", fa[1])
  } else {
    r <- am |> distinct(position, aa_ref)
    stopifnot(!anyDuplicated(r$position))
    s <- rep(NA_character_, max(r$position))
    s[r$position] <- r$aa_ref
    if (is.na(s[1])) s[1] <- "M"     # start codon is not scored by alphamissense
    message("reference rebuilt from the alphamissense table; gaps: ", sum(is.na(s)))
  }
  if (length(s) != 422) warning("reference length is ", length(s), ", expected 422")
  s
}

message("setup ok")