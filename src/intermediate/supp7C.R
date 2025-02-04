library("dplyr")
library("readr")
library("tibble")
library("tidyr")
library("DESeq2")
library("glue")

dir <- "~/huo2025/data/rnaseq/"
cts_file <- "processed_reads_for_deseq2.csv"
col_file <- "input_experiment.csv"
fig_dir <- "~/huo2025/fig/"
# Load
cts <- read_csv(file.path(dir, cts_file)) %>% column_to_rownames("Geneid")
coldata <- read_csv(file.path(dir, col_file)) 
# Wrangle
coldata <- coldata %>% column_to_rownames("SampleID")
tp.to.midlog <- Vectorize(function(x) {if (x == "0h") {"mid-log"} else {paste("mid-log", x, sep="\n+ ")}})
fru.frupal.expand <- Vectorize(function (x) {if (x == "Fru") {"0.1% Fru"} else {"0.1% Fru, \n0.9% Pal"}})
add.delta <- Vectorize(function (x) {if (x == "gal80") {"gal80Δ"} else {x}})
coldata <- coldata %>% mutate(Timepoint = tp.to.midlog(Timepoint), Condition = fru.frupal.expand(Condition), Genotype = add.delta(Genotype))
sample_names <- rownames(coldata)
cts <- cts %>% select(any_of(sample_names))
cts <- as.matrix(cts)
coldata$Genotype <- factor(coldata$Genotype, levels = c("WT", "gal80Δ"))
coldata$Condition <- factor(coldata$Condition, levels = c("0.1% Fru", "0.1% Fru, \n0.9% Pal"))
coldata$Timepoint <- factor(coldata$Timepoint, levels = c("mid-log", "mid-log\n+ 10h", "mid-log\n+ 16h"))
coldata$BioRep <- factor(coldata$BioRep, levels= c("A", "B", "C"))
coldata
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Genotype + Condition + Timepoint)
design <- c("Genotype", "Condition", "Timepoint")
vsd <- vst(dds, blind=TRUE)

ppi <- 300
png(file.path(fig_dir, "supp7c.png"), width = 3.5*ppi, height = 4*ppi, res = ppi)
pcaData <- plotPCA(vsd, intgroup=design, returnData=TRUE, ntop=dim(vsd)[[1]])
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=Genotype, shape=Timepoint, alpha=Condition)) +
  geom_point(size=5) +
  scale_alpha_manual(values=c(0.25, 0.75)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(legend.position = "bottom", legend.title.position = "top", legend.text = element_text(size = 8), legend.key.size = unit(0.75, "lines")) + 
  guides(colour = guide_legend(nrow = 3), shape = guide_legend(nrow = 3), alpha = guide_legend(nrow = 2))
dev.off()
