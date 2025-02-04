library(tidyqpcr)
library(dplyr)
library(repr)
library(ggplot2)

dir <- "~/huo2025/data/supp7/"
fig_dir <- "~/huo2025/fig/"

gene_name_levels <- c("GAL1", "GAL2", "GAL4", "GAL7", "GAL80", "ACT1", "ALG9", "PUS7")
gene_name_values <- rep(gene_name_levels, times = 2)


strain_values_1 <- c(rep(c("gal80", "FY4"), each=5),
                   c("FY4"))
strain_values <- c(rep(strain_values_1, times=2), rep(c(""),2))
RT_status <- c(rep(c("+RT", "-RT"), each=11), rep(c("NT"), 2))
condition_status_1 <- c("Gal",
                      rep(c("Fru", "Gal"), times=5))
condition_status <- c(rep(condition_status_1, times=2), rep(c(""), 2))
biol_rep_values_1 <- c(rep(c("A", "B", "B" ,"C", "C")),
                     rep(c("A", "B", "C"), each=2))
biol_rep_values <- c(rep(biol_rep_values_1, times=2), rep(c(""), 2))


rowkey <- tibble(
  well_row = LETTERS[1:16],
  target_id = gene_name_values
)

colkey <- tibble(
  well_col = c(1:24),
  strain = strain_values,
  prep_type = RT_status,
  condition = condition_status,
  biol_rep = biol_rep_values,
) %>% mutate(
        sample_id = paste(strain, biol_rep, condition, sep = "_")
      )

plate1plan <-
  label_plate_rowcol(
    create_blank_plate(),
    rowkey,
    colkey
  )

setwd(dir)
plate_cq_data <- read_lightcycler_1colour_cq(
  "20220217_Yu_qPCR_GAL_in_fru_and_gal_Cq.txt",
) %>% right_join(plate1plan)

plate_norm <- plate_cq_data[plate_cq_data$prep_type == "+RT",] %>%
  calculate_deltacq_bysampleid(ref_target_ids = c("PUS7", "ALG9", "ACT1"), norm_function = median)

plate_deltanorm <- plate_norm %>%
  calculate_deltadeltacq_bytargetid(ref_sample_ids = c("FY4_A_Fru", "FY4_B_Fru", "FY4_C_Fru"))

plate_deltanorm_med <- plate_deltanorm %>%
    group_by(sample_id, strain, condition, target_id) %>%
  summarize(
    deltadelta_cq  = median(deltadelta_cq, na.rm = TRUE),
    fold_change    = median(fold_change,   na.rm = TRUE)
  )

ppi <- 300
png(file.path(fig_dir, "supp7a.png"), width = 3*ppi, height = 4*ppi, res = ppi)
ggplot(data = filter(plate_deltanorm_med, !(target_id %in% c("PUS7", "ALG9", "ACT1")))) +
  geom_point(aes(x = condition, y = deltadelta_cq, shape = strain, colour = strain), size=5, alpha=0.7) +
  labs(
    y = "ΔΔCq (log2 fold change)\n relative to FY4 in fructose"
  ) +
  facet_wrap(~target_id,ncol=6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
dev.off()
