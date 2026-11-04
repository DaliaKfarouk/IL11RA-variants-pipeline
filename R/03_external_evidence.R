# external evidence read straight from the exports (nothing typed by hand)
source("R/00_setup.R")
dat  <- readRDS("results/m1_02_predicted.rds")
am   <- read_am()
f_gn <- locate("^gnomad.*\\.csv$")
f_cv <- locate("^clinvar.*\\.(tsv|txt)$")
f_am <- locate("alphamissense.*\\.(zip|tsv)$")

# version and date trail
tibble(resource = c("gnomad", "clinvar", "alphamissense"),
       file = basename(c(f_gn, f_cv, f_am)),
       md5  = unname(tools::md5sum(c(f_gn, f_cv, f_am))),
       file_date = as.character(as.Date(file.mtime(c(f_gn, f_cv, f_am))))) |>
  write_csv("results/m1_external_sources.csv")
gn_version <- str_extract(basename(f_gn), "v[0-9]+(_[0-9]+)*") |> str_replace_all("_", ".")

# ---- gnomad: join on chrom / pos / ref / alt ----
gn <- read_csv(f_gn, col_types = cols(.default = "c"), na = c("", "NA"), show_col_types = FALSE)
check_cols(gn, c("Chromosome", "Position", "Reference", "Alternate", "Allele Count",
                 "Allele Number", "Allele Frequency", "Homozygote Count"), "gnomad")
gn_k <- gn |>
  transmute(chrom_bare = str_remove(Chromosome, "^chr"), pos = as.integer(Position),
            ref = Reference, alt = Alternate,
            gnomad_ac = as.numeric(`Allele Count`), gnomad_an = as.numeric(`Allele Number`),
            gnomad_af = as.numeric(`Allele Frequency`),
            gnomad_hom = as.numeric(`Homozygote Count`)) |>
  distinct(chrom_bare, pos, .keep_all = TRUE)

dat <- dat |>
  mutate(chrom_bare = str_remove(chrom_ucsc, "^chr")) |>
  left_join(gn_k, by = c("chrom_bare", "pos", "ref", "alt")) |>
  mutate(gnomad_version = gn_version, gnomad_matched = !is.na(gnomad_ac))
message(sum(dat$gnomad_matched), "/", nrow(dat), " variants found in the gnomad export")

# ---- clinvar: join on the cDNA change inside the Name column ----
cv <- read_tsv(f_cv, col_types = cols(.default = "c"), na = c("", "NA"), show_col_types = FALSE)
check_cols(cv, c("Name", "Accession", "Germline classification", "Germline review status"), "clinvar")
cv_k <- cv |>
  transmute(hgvs_c = str_extract(Name, "c\\.[^\\s)]+"), clinvar_id = Accession,
            clinvar_class = `Germline classification`,
            clinvar_review_status = `Germline review status`) |>
  filter(!is.na(hgvs_c)) |> distinct(hgvs_c, .keep_all = TRUE)
dat <- left_join(dat, cv_k, by = "hgvs_c")
message(sum(!is.na(dat$clinvar_id)), "/", nrow(dat), " variants matched in clinvar")
if (all(is.na(dat$clinvar_id))) print(head(cv$Name, 5))

# ---- alphamissense: missense only, key = protein change ----
am_k <- am |> transmute(am_key = protein_variant, am_score, am_class) |>
  distinct(am_key, .keep_all = TRUE)
dat <- dat |>
  mutate(am_key = if_else(var_class == "missense", coalesce(predicted_hgvs_p, hgvs_p_n), NA_character_)) |>
  left_join(am_k, by = "am_key") |>
  mutate(am_status = case_when(!is.na(am_class) ~ "scored",
                               var_class == "missense" ~ "missense_not_found",
                               TRUE ~ "not_scored"))
print(count(dat, am_status))

saveRDS(dat, "results/m1_03_evidence.rds")