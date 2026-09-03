library(data.table)
library(dplyr)
library(BuenColors)

aou <- fread("../../anello-covariates/biobanks/all-of-us/data/non_zero_entries_anello.tsv")
ukb <- fread("../../anello-covariates/biobanks/uk-biobank/crude/UKB_non_zero_entries_kallisto.tsv")
meta <- fread("../../anello-covariates/color-scheme/final_simple_anello_metadata.csv")

aou_filtered <- aou %>% filter(est_counts >= 1) %>%
  group_by(target_id) %>%
  summarize(n_donors = n(), total_reads = sum(est_counts)) %>%
  arrange(desc(n_donors))

ukb_filtered <- ukb %>% filter(est_counts >= 1) %>%
  group_by(target_id) %>%
  summarize(n_donors = n(), total_reads = sum(est_counts)) %>%
  arrange(desc(n_donors))

# 1. Combine UKB and AoU into one long format
combined_summary <- bind_rows(ukb_filtered, aou_filtered) %>%
  # 2. Group by the virus ID and sum the metrics
  group_by(target_id) %>%
  summarise(
    n_donors_total = sum(n_donors, na.rm = TRUE),
    total_reads_sum = sum(total_reads, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 3. Bring in the metadata from the 'meta' df
  # This matches Accession (meta) to target_id (the viral IDs)
  left_join(meta, by = c("target_id" = "Accession"))

head(combined_summary)

#taking the top and bottom prevalence anello strains by genus for ORF1/2 expression and RNAseq experiment
target_genera <- c("Alphatorquevirus", "Betatorquevirus", "Gammatorquevirus")

extreme_anellos <- combined_summary %>%
  # 1. Convert to a standard base-R dataframe to strip 'Rle' or 'S4' classes
  as.data.frame() %>%
  # 2. Basic filter
  filter(orf1_genus %in% target_genera, n_donors_total >= 10) %>%
  # 3. Use group_by and a simple rank filter
  group_by(orf1_genus) %>%
  arrange(dplyr::desc(n_donors_total), .by_group = TRUE) %>%
  # 4. Take the first and last row of each group
  mutate(selected = row_number() == 1 | row_number() == n()) %>%
  mutate(rank = 1:n()) %>%
  ungroup() %>%
  # 5. Final sort
  arrange(orf1_genus, dplyr::desc(n_donors_total))

p1 <- ggplot(extreme_anellos, aes(x = rank, y = n_donors_total, color = selected)) +
  geom_point(size = 2) + scale_y_log10() + facet_wrap(~orf1_genus) + scale_x_log10() +
  pretty_plot(fontsize = 7) + 
  labs(x = "Rank", y = "n Donors (AoU+UKB)")+
  theme(legend.position = "none") +
  scale_color_manual(values = c("lightgrey", "firebrick"))

cowplot::ggsave2(p1, file = "../plots/which_strains.pdf", width = 3.2, height = 1.8)
