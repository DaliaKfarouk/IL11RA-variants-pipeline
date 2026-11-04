# predictCoding (VariantAnnotation vignette workflow), concordance, reference residue check
source("R/00_setup.R")

dat     <- readRDS("results/m1_01_prepared.rds")
am      <- read_am()
ref_seq <- uniprot_ref(am)
txdb    <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome  <- BSgenome.Hsapiens.UCSC.hg38

# predictCoding is built for substitutions: snvs only, other classes keep the SR class
snv <- dat |> filter(is_snv, !is.na(pos), !is.na(chrom_ucsc)) |> mutate(qid = row_number())
stopifnot("no SNVs with coordinates" = nrow(snv) > 0)

gr  <- GRanges(seqnames = snv$chrom_ucsc, ranges = IRanges(start = snv$pos, width = 1))
cod <- predictCoding(query = gr, subject = txdb, seqSource = genome,
                     varAllele = DNAStringSet(snv$alt))   # alt on the + strand

mc <- mcols(cod)
cand <- tibble(
  qid      = as.integer(mc$QUERYID),
  txid     = as.character(mc$TXID),
  conseq   = as.character(mc$CONSEQUENCE),
  prot_pos = vapply(as.list(mc$PROTEINLOC),
                    \(x) if (length(x)) as.integer(x[[1]]) else NA_integer_, integer(1)),
  ref_aa   = na_if(unname(as.character(mc$REFAA)), ""),
  var_aa   = na_if(unname(as.character(mc$VARAA)), ""))

cds <- cdsBy(txdb, by = "tx")
cds_len <- tibble(txid = names(cds), cds_len = as.integer(sum(width(cds))))
tx_name <- AnnotationDbi::select(txdb, keys = unique(cand$txid), columns = "TXNAME",
                                 keytype = "TXID") |>
  as_tibble() |> transmute(txid = as.character(TXID), tx_name = TXNAME)

# isoform rule: one row per variant. prefer the transcript whose reference residue matches
# uniprot Q14626, then the longest CDS. n_tx and isoforms_agree keep the isoform differences
cand <- cand |>
  left_join(cds_len, by = "txid") |> left_join(tx_name, by = "txid") |>
  mutate(ref_ok = ref_aa == ref_seq[prot_pos])

pick <- cand |>
  group_by(qid) |>
  mutate(n_tx = n(), isoforms_agree = n_distinct(paste(prot_pos, ref_aa, var_aa)) == 1) |>
  arrange(desc(ref_ok), desc(cds_len), txid, .by_group = TRUE) |>
  slice(1) |> ungroup()

aa3 <- function(x) unname(AA3[x])
pred <- snv |>
  left_join(pick, by = "qid") |>
  mutate(predicted_hgvs_p = case_when(
    conseq == "synonymous" ~ paste0("p.", aa3(ref_aa), prot_pos, "="),
    !is.na(var_aa)         ~ paste0("p.", aa3(ref_aa), prot_pos, aa3(var_aa)),
    TRUE                   ~ NA_character_)) |>
  select(variant_id, predicted_consequence = conseq, predicted_hgvs_p,
         predicted_tx = tx_name, n_tx, isoforms_agree)

dat <- dat |>
  left_join(pred, by = "variant_id") |>
  mutate(hgvs_p_n = norm_p(hgvs_p), mane_p_n = norm_p(mane_hgvs_p),
         protein_match = predicted_hgvs_p == hgvs_p_n,       # vs the SR
         mane_match    = predicted_hgvs_p == mane_p_n)       # vs VariantValidator

# reference residue check from the published protein name (all rows with a p. name)
m <- str_match(dat$hgvs_p_n, "^p\\.([A-Z][a-z]{2})([0-9]+)")
dat$prot_pos_published   <- as.integer(m[, 3])
dat$ref_aa_published     <- unname(AA1[m[, 2]])
dat$reference_residue_ok <- dat$ref_aa_published == ref_seq[dat$prot_pos_published]

chk <- dat |> filter(!is.na(protein_match))
message(sprintf("protein change concordant: %d/%d (%.1f%%)",
                sum(chk$protein_match), nrow(chk), 100 * mean(chk$protein_match)))
dat |> filter(protein_match %in% FALSE) |>
  select(variant_id, hgvs_c, hgvs_p, mane_hgvs_p, predicted_hgvs_p, predicted_tx, n_tx, isoforms_agree) |>
  write_csv("results/m1_concordance_mismatches.csv")
bad_ref <- dat |> filter(reference_residue_ok %in% FALSE) |>
  select(variant_id, hgvs_p, prot_pos_published, ref_aa_published)
if (nrow(bad_ref)) { print(bad_ref); message("reference residue mismatches above: offset or isoform?") }

saveRDS(dat, "results/m1_02_predicted.rds")