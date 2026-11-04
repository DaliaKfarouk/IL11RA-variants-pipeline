# module 2: domain table, position table, lollipop, domain counts
source("R/00_setup.R")
dat <- readRDS("results/m1_04_filters.rds")

dom <- read_csv("data/domains_IL11RA.csv", col_types = cols(.default = "c"),
                na = c("", "NA"), show_col_types = FALSE) |>
  mutate(start = as.integer(start), end = as.integer(end))
stopifnot("fill start/end in data/domains_IL11RA.csv from UniProt (delete unused rows)" =
            !anyNA(dom$start) && !anyNA(dom$end))

pts <- dat |>
  filter(!is.na(prot_pos_published)) |>
  mutate(func_class = case_when(
    is.na(functional_evidence) | str_detect(functional_evidence, regex("^none|not (tested|assessed)", TRUE)) ~ "none",
    str_detect(functional_evidence, regex("retention|traffick|surface", TRUE)) ~ "impaired trafficking",
    str_detect(functional_evidence, regex("pstat3|signal", TRUE))              ~ "reduced signalling",
    TRUE ~ "other"),
    low_quality = is.na(quality_rating) | !str_detect(quality_rating, regex("^low risk", TRUE)))
message(nrow(dat) - nrow(pts), " variants without a protein position are not on the map")

# one point per position, y = independent families (a family reported twice counts once)
pos_tbl <- pts |>
  group_by(prot_pos_published, var_class, func_class) |>
  summarise(n_families = n_distinct(family_id, na.rm = TRUE),
            any_low_quality = any(low_quality),
            label = paste(unique(hgvs_p_n), collapse = "/"), .groups = "drop")
write_csv(pos_tbl, "results/m2_position_table.csv")

dom_counts <- dom |>
  mutate(length_aa = end - start + 1) |>
  rowwise() |>
  mutate(n_variants = sum(between(pts$prot_pos_published, start, end)),
         n_families = n_distinct(pts$family_id[between(pts$prot_pos_published, start, end)], na.rm = TRUE),
         variants_per_100aa = round(100 * n_variants / length_aa, 1)) |>
  ungroup()
write_csv(dom_counts, "results/m2_domain_counts.csv")
print(dom_counts)

y0 <- -0.4
p <- ggplot() +
  geom_rect(data = dom, aes(xmin = start, xmax = end, ymin = y0 - 0.3, ymax = y0 + 0.3),
            fill = "grey88", colour = "grey40") +
  geom_text(data = dom, aes(x = (start + end) / 2, y = y0, label = feature), size = 3) +
  geom_segment(data = pos_tbl, aes(x = prot_pos_published, xend = prot_pos_published,
                                   y = 0, yend = n_families), colour = "grey50") +
  geom_point(data = pos_tbl, aes(prot_pos_published, n_families, colour = var_class,
                                 shape = func_class), size = 3.5) +
  geom_text_repel(data = pos_tbl,
                  aes(prot_pos_published, n_families,
                      label = paste0(label, if_else(any_low_quality, "*", ""))), size = 3) +
  scale_y_continuous(breaks = function(l) seq(0, ceiling(l[2]), by = 1)) +
  labs(x = "residue (UniProt Q14626 numbering)", y = "independent families",
       colour = "variant class", shape = "functional evidence",
       caption = "* report not rated low risk of bias. source: SR variant table; domains: UniProt") +
  theme_minimal()
ggsave("figures/m2_lollipop_IL11RA.png", p, width = 10, height = 4.5, dpi = 300)