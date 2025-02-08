library(tidyqpcr)
library(dplyr)
library(repr)
library(ggplot2)
library(rstudioapi)

script_path <- rstudioapi::getSourceEditorContext()$path

dir <- dirname(dirname(dirname(script_path)))
data_dir <- file.path(dir, "data/supp6")
fig_dir <- file.path(dir, "fig")

# Names of target genes
gene_name_values <- c(rep(c("GAL1", "GAL2", "GAL4", "GAL7", "GAL80", "MAL11"), times=2), rep(c("MAL13", "ZNF1", "IMA5", "ACT1", "ALG9", "PUS7"),times=2))
strain_values_1 <- rep(c("FY4", "gal80", "GAL2-OE"), each=2)
strain_values <- c(rep(strain_values_1, times=2), c(""))
RT_status <- c(rep(c("+RT", "-RT"), each=6), c("NT"))
condition_status <- c(rep("Fru",times=12), c(""))
biol_rep_values <- c(rep(c("A", "B"), times=6), c(""))

rowkey <- tibble(
  well_row = LETTERS[1:13],
  strain = strain_values,
  prep_type = RT_status,
  condition = condition_status,
  biol_rep = biol_rep_values,
) %>% mutate(
        sample_id = paste(strain, biol_rep, condition, sep = "_")
      )

colkey <- tibble(
  well_col = c(1:24),
  target_id = gene_name_values
  )

plate1plan <-
  label_plate_rowcol(
    create_blank_plate(),
    rowkey,
    colkey
  )

setwd(data_dir)
plate_cq_data <- read_lightcycler_1colour_cq(
  "20220615_Yu_qPCR_GAL2_OE_Cq.txt",
) %>% right_join(plate1plan)


# Removed all -RT and NT controls before normalising!
plate_norm <- plate_cq_data[plate_cq_data$prep_type == "+RT",] %>%
  calculate_deltacq_bysampleid(ref_target_ids = c("PUS7", "ALG9", "ACT1"), norm_function = mean)

plate_deltanorm <- plate_norm %>%
  calculate_deltadeltacq_bytargetid(ref_sample_ids = c("FY4_A_Fru", "FY4_B_Fru"), norm_function = mean)

plate_deltanorm_med <- plate_deltanorm %>%
    group_by(sample_id, strain, condition, target_id) %>%
  summarize(
    deltadelta_cq  = mean(deltadelta_cq, na.rm = TRUE),
    fold_change    = mean(fold_change,   na.rm = TRUE)
  )



ppi <- 300
png(file.path(fig_dir, "supp6c.png"), width = 3.5*ppi, height = 4*ppi, res = ppi)
ggplot(data = filter(plate_deltanorm_med, target_id %in% c("GAL1", "GAL2", "GAL80"))) +
  geom_point(aes(x = condition, y = deltadelta_cq, colour = strain), size=5, alpha=0.7) +
  labs(
    y = "ΔΔCq (log2 fold change)\n relative to FY4 in fructose"
  ) +
  facet_wrap(~target_id,ncol=6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), legend.position = "bottom")
dev.off()

