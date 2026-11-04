# freeze input, join coordinates, define variant classes
source("R/00_setup.R")

in_path <- locate("^variants_input\\.csv$")
co_path <- locate("^variants_coords\\.csv$")
rd <- function(f) read_csv(f, col_types = cols(.default = "c"), na = c("", "NA"), show_col_types = FALSE)
variants <- rd(in_path)
coords   <- rd(co_path)

check_cols(variants, c("variant_id", "gene", "hgvs_c", "hgvs_p", "variant_type", "family_id",
                       "pmid", "quality_rating", "functional_evidence"), "variants_input")
cols_c <- c("variant_id", "genome_build", "chrom", "pos", "ref", "alt",
            "mane_hgvs_c", "mane_hgvs_p", "validation_status")
check_cols(coords, cols_c, "variants_coords")
stopifnot(!anyDuplicated(variants$variant_id), !anyDuplicated(coords$variant_id))
stopifnot("coordinates must be GRCh38" = all(str_detect(na.omit(coords$genome_build), "38")))

# record what was frozen
tibble(file = c(in_path, co_path), md5 = unname(tools::md5sum(c(in_path, co_path))),
       date = as.character(Sys.Date())) |>
  write_csv("results/m1_input_freeze.csv")

# other genes are kept apart, rows without any name are not analysable
variants |> filter(gene != "IL11RA" | is.na(gene)) |> write_csv("results/m1_other_genes.csv")

dat <- variants |>
  filter(gene == "IL11RA", !is.na(hgvs_c) | !is.na(hgvs_p)) |>
  left_join(select(coords, all_of(cols_c)), by = "variant_id") |>
  mutate(pos = as.integer(pos),
         chrom_ucsc = if_else(is.na(chrom), NA_character_, paste0("chr", str_remove(chrom, "^chr"))),
         is_snv    = coalesce(str_detect(ref, "^[ACGT]$") & str_detect(alt, "^[ACGT]$"), FALSE),
         validated = coalesce(validation_status == "validated", FALSE),
         canonical_splice = coalesce(str_detect(hgvs_c, "[0-9][+-][12](?![0-9])"), FALSE),
         var_class = case_when(
           str_detect(variant_type, regex("frameshift", TRUE))                   ~ "frameshift",
           str_detect(variant_type, regex("nonsense|stop", TRUE))                ~ "nonsense",
           str_detect(variant_type, regex("splice", TRUE)) & canonical_splice    ~ "canonical_splice",
           str_detect(variant_type, regex("splice", TRUE))                       ~ "splice_region",
           str_detect(variant_type, regex("in-?frame|indel|duplication|insertion|deletion", TRUE)) ~ "inframe_indel",
           str_detect(variant_type, regex("missense", TRUE))                     ~ "missense",
           str_detect(variant_type, regex("synonymous", TRUE))                   ~ "synonymous",
           TRUE ~ "other"))

# check the class mapping by eye before going on
dat |> count(variant_type, var_class) |> print(n = Inf)
message(nrow(dat), " IL11RA rows | ", sum(!dat$validated), " unvalidated (kept) | ", sum(dat$is_snv), " SNVs")

saveRDS(dat, "results/m1_01_prepared.rds")