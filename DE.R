library(DESeq2)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tidyr)

# reading STAR gene counts
read.csv2('ins2ReadsPerGene.out.tab', header=FALSE, sep='\t', 
                  stringsAsFactors = FALSE)  %>% 
  select(gene = 1, ins2 = 2)-> ins2
read.csv2('a1s7-mix-filtReadsPerGene.out.tab', header=FALSE, sep='\t',
                  stringsAsFactors = FALSE)  %>% 
  select(gene = 1, a1s7 = 2) -> a1s7
read.csv2('orts1-mix-filtReadsPerGene.out.tab', header=FALSE, sep='\t',
                   stringsAsFactors = FALSE) %>% 
  select(gene = 1, orts1 = 2) -> orts1

# reading metadata file
MetaData <- read.csv2('MetaData', header=TRUE, sep='\t')

# merging counts in one dataset and removing L. longissimus evidence
merge(a1s7[0:8724,], orts1[0:8724,], by='gene') %>% 
  merge(ins2[0:8724,], by='gene') %>%
  tibble::column_to_rownames('gene') -> CountMatrix

# writing count matrix to file
write.csv2(CountMatrix, file="CountMatrix.csv")  

# creating DESeq2 object 
dds <- DESeqDataSetFromMatrix(countData = CountMatrix,
                              colData = MetaData,
                              design = ~ condition)
dds <- DESeq(dds)

# plotting PCA to check sample variance
vst <- varianceStabilizingTransformation(dds)
plotPCA(vst)

# calculating DE results and writing them to file
res <- results(dds)
summary(res)
res <- res[order(res$padj),]
write.csv2(res, file="res-all.csv")

########## EdgeR expression threshold ##########

# EdgeR metadata
MetaData <- c("free", "free", "parasite")

# creating edgeR DGE object from existing CountMatrix and EdgeR Metadata
dge <- DGEList(counts = CountMatrix, group = factor(MetaData))

# filtering rows with 0 across all samples and calculating normalized factors
nonzero <- rowSums(dge$counts) > 0
dge %>% .[nonzero,] %>% calcNormFactors -> dge

#nrow(dge)
#[1] 8570

# average logCPM density plot for free-living samples
a <- aveLogCPM(dge[,1:2])
dens <- density(a)
plot(dens, las=1, main="LogCPM distribution in free-living samples")

# manually identifying density threshold (the point is between two density peaks)
tr <- locator(n=1)
legend(x = tr[1], y= tr[2], legend = round(tr$x, 3), bty="n")
avelogcpm.presence.threshold <- tr$x
# avelogcpm.presence.threshold = -1.829106
# CPM = 0.281 

# plotting average logCPM density distribution with expression threshold
ggplot(data.frame(logCPM=a)) +
  aes(x=logCPM) +
  geom_histogram(aes(y=100*(..count..)/sum(..count..)), binwidth=0.25, boundary=0) +
  geom_vline(xintercept=avelogcpm.presence.threshold, color="red", linetype="dashed") +
  xlab("Average logCPM") + ylab("Percent of genes in bin") +
  coord_cartesian(xlim=quantile(a, c(0, 0.995)), ylim=c(0,4)) +
  labs(title="Average gene LogCPM distribution",
       subtitle="for genes with at least 1 read") +
  theme(plot.caption = element_text(hjust = 0))

# selecting genes unexpressed in free-living samples (below expression threshold)
dge[,1:2][a < avelogcpm.presence.threshold,] -> dgefree
edger.unexpr <- as.data.frame(dgefree$counts)
write.csv2(edger.unexpr, file="dge-unexpr-in-free-new1.csv")

########## back to DESeq2 ########## 

# keeping only DE genes unexpressed in free-living samples 
# with p.adj <= 0.05 and log2FC > 2
as.data.frame(res) %>% 
  filter((row.names(res) %in% c(row.names(edger.unexpr)))) %>%
  filter(padj < 0.05  & log2FoldChange > 2)-> res.unexpr.deed
write.csv2(res.unexpr.deed, file="res.unexpr.deed.csv")