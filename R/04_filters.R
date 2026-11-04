# F1 / F2 / F3, sensitivity with exact binomial CIs, filter matrix figure
source("R/00_setup.R")
dat <- readRDS("results/m1_03_evidence.rds")

# thresholds: same numbers go in docs/analysis_plan.md
F1_AF <- 0.01
F2_AF <- 0.001
F3_AF <- 0.001
NA_AF_PASSES <- TRUE    # absent from gnomad counts as rare (ACMG PM2 logic)

af_ok <- function(af, thr) if_else(is.na(af), NA_AF_PASSES, af < thr)
why   <- function(pass, af_pass, cls_pass)
  case_when(pass ~ "pass", !af_pass ~ "frequency", !cls_pass ~ "variant class / score", TRUE ~ "other")

dat <- dat |>
  mutate(lof    = var_class %in% c("nonsense", "frameshift", "canonical_splice"),
         mis_lp = var_class == "missense" & am_class %in% "likely_pathogenic",
         cls1 = var_class != "synonymous",
         cls2 = lof | mis_lp,
         cls3 = cls2 | var_class %in% c("splice_region", "inframe_indel", "missense"),
         passes_F1 = af_ok(gnomad_af, F1_AF) & cls1,
         passes_F2 = af_ok(gnomad_af, F2_AF) & cls2,
         passes_F3 = af_ok(gnomad_af, F3_AF) & cls3,
         why_F1 = why(passes_F1, af_ok(gnomad_af, F1_AF), cls1),
         why_F2 = why(passes_F2, af_ok(gnomad_af, F2_AF), cls2),
         why_F3 = why(passes_F3, af_ok(gnomad_af, F3_AF), cls3))

filt <- dat |>
  select(variant_id, hgvs_p, hgvs_c, var_class, starts_with("passes_")) |>
  pivot_longer(starts_with("passes_"), names_to = "filter", values_to = "pass") |>
  mutate(filter = str_remove(filter, "passes_"))

# sensitivity = share of known disease variants that pass (all SR variants are disease-reported)
sens <- filt |>
  group_by(filter) |>
  summarise(n = n(), passed = sum(pass), .groups = "drop") |>
  mutate(sensitivity = passed / n,
         ci_low  = map2_dbl(passed, n, \(x, k) binom.test(x, k)$conf.int[1]),
         ci_high = map2_dbl(passed, n, \(x, k) binom.test(x, k)$conf.int[2]))
print(sens)
write_csv(sens, "results/m1_filter_sensitivity.csv")

dat |> group_by(var_class) |>
  summarise(n = n(), F1 = sum(passes_F1), F2 = sum(passes_F2), F3 = sum(passes_F3)) |>
  write_csv("results/m1_filter_by_class.csv")

dat |> filter(!(passes_F1 & passes_F2 & passes_F3)) |>
  select(variant_id, hgvs_c, hgvs_p, var_class, gnomad_af, am_class, starts_with("why_")) |>
  write_csv("results/m1_filter_missed.csv")

filt |> pivot_wider(names_from = filter, values_from = pass, names_prefix = "passes_") |>
  write_csv("results/m1_filter_matrix.csv")

# source concordance: clinvar vs alphamissense, described not resolved
dat |> count(clinvar_class, am_class, name = "n") |> write_csv("results/m1_source_concordance.csv")

p <- filt |>
  mutate(label = paste(variant_id, coalesce(hgvs_p, hgvs_c)),
         pass  = factor(pass, c(TRUE, FALSE), c("pass", "fail")),
         filter = toupper(filter)) |>
  ggplot(aes(filter, fct_reorder(label, as.integer(factor(var_class))), fill = pass)) +
  geom_tile(colour = "white") +
  scale_fill_manual(values = c(pass = "#1b9e77", fail = "#d95f02")) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal()
ggsave("figures/m1_filter_matrix.png", p, width = 5, height = 8, dpi = 300)

dat |> select(-any_of(c("lof", "mis_lp", "cls1", "cls2", "cls3"))) |>
  write_csv("results/m1_variants_annotated.csv")
saveRDS(dat, "results/m1_04_filters.rds")