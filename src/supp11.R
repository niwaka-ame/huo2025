library("dplyr")
library("readr")
library("tibble")
library("tidyr")
library("pheatmap")
library("glue")

dir <- "~/huo2025/data/rnaseq/"
fig_dir <- "~/huo2025/fig/"

setwd(dir)
for (i in c(0, 10 ,16)) {
  fname <- glue("group_gal80.Fru.{i}h_vs_WT.Fru.{i}h.csv")
  if (i == 0){
    tb <- read_csv(fname) %>% filter(padj<=0.05, abs(log2FoldChange)>=0.5)
    tb$Timepoint <- glue("mid-log") 
    tb$Timepoint <- tb$Timepoint %>% factor()
    genes <- unique(tb$...1)
  } else {
    tb_temp <- read_csv(fname) %>% filter(padj<=0.05, abs(log2FoldChange)>=0.5, ...1 %in% genes)
    tb_temp$Timepoint <- glue("mid-log\n+ {i} h")
    tb_temp$Timepoint <- tb_temp$Timepoint %>% factor()
    genes <- unique(tb_temp$...1) # update genes
    tb <- full_join(tb, tb_temp)
  }
}
tb <- tb %>% filter(...1 %in% genes) %>% select(...1, log2FoldChange, Timepoint)

dfgene <- read_table(file = "yeast.txt", skip=58, col_names=FALSE) %>% rename(gene = X1, systematic = X2) %>% select(gene, systematic)

split.genename <- function(g, sep) unlist(strsplit(g, sep))[1]
split.genename <- Vectorize(split.genename)
tb <- tb %>% mutate(systematic = split.genename(...1, "_")) %>% left_join(dfgene, by = "systematic") %>% select(gene, log2FoldChange, Timepoint) %>% mutate(gene = split.genename(gene, ";"))

hm.mat <- tb %>% pivot_wider(names_from = Timepoint, values_from = log2FoldChange) %>% column_to_rownames(var="gene") %>% as.matrix()

# breaks <- seq(-2, 2, length.out = 101)
# pheatmap(hm.mat, fontsize_row=4, breaks= breaks, cellheight = 4, cellwidth = 20, cutree_rows = 5)

# For presentation
breaks <- seq(-2, 2, length.out = 101)
ppi <- 600
png(file.path(fig_dir, "supp11.png"), height = 6*ppi, width = 3*ppi, res = ppi)
pheatmap(hm.mat, fontsize_row=4, breaks= breaks, cellheight = 4, cellwidth = 20, cutree_rows = 5, angle_col=0, fontsize_col = 5)
dev.off()
