###Code Section 1. Data Acquisition, exploration, filtering and quality

library("airway")#Dataset (Himes et al., 2014)

#The primary goal of this workflow was to generate a start to end pipeline for gene expression analysis similar to what the RNA-seq workflow vignette by Love et al., 2019; starting with the estimation of transcript-level abundance from raw RNA seq read data with a program as Salmon. Unfortunately, Salmon was tested on the command line interface of Windows Subsystem for Linux in my case but the analysis is memory demanding and computationally intensive, so the present workflow resumes on the output of Salmon which is provided by the "Airway" package. 

#Demonstrative section on how to import and process files directly outputed by Salmon. It is demonstrative as it only takes two files of the whole dataset.

##Loading Airway dataset
a.dir <- system.file("extdata", package="airway", mustWork=TRUE) #printing local directory where Airway files are
list.files(file.path(a.dir, "quants")) #List all read data files from quants dir (only 2 files)

##Table (csv) linking sample details and location of fastq files (needed for tximeta)
#This file needs to be manually generated
csv.file <- file.path(a.dir, "sample_table.csv") 
col.data <- read.csv(csv.file, row.names=1, stringsAsFactors=FALSE) #reading table (csv)
head(col.data, n = 3L) #Print first 3 rows of table

#Modifying table to only work with two read data files as just "SRR1039508" "SRR1039509" are included.
col.data <- col.data[1:2,] 
col.data$names <- col.data$Run
col.data$files <- file.path(a.dir, "quants", col.data$names, "quant.sf.gz") #Copying location of files
file.exists(col.data$files)

#Tximeta-importing data into R
library("tximeta") #Loves et al., 2020
#Generating transcript-level quantifications
t.data <- tximeta(col.data)
#From the output, it is important that the program matches the transcriptome used, in this case: GENCODE - Homo sapiens - release 29

#Datacheck
dim(t.data)
head(rownames(t.data))

#Summarizing transcript-level quantitifications to gene levels
g.data <- summarizeToGene(t.data)
#Datachecks
dim(g.data)
head(rownames(g.data))

#At this point ends the demonstrative section
#Now we load the complete summarized data for 8 samples (already in package "Airway" in the "gse" file)
data(gse)
dim(gse) 
head(rownames(gse))
assayNames(gse)

gse$donor
gse$condition

#Renaming variables
#Cell as cell line and dex as treatment condition
gse$cell <- gse$donor
gse$dex <- gse$condition

#Changing name of levels
#Untreated as untrt and Dexamethasone as trt
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")
levels(gse$dex)

#Checking millions of fragments that could be mapped by Salmon to the genes, second arguments indicates number of decimal points. 
round(colSums(assay(gse)) / 1e6, 1 )

#DESeq2 Data Load
#Using counts and average transcript lengths from tximeta
#Design formula for differential gene expression
# ~condition would be the simplest one, in this case we will use ~ cell + dex
#Testing the effect of dexamethasone controlling for the effect of differential cell line
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ cell + dex)
class(dds)
dim(dds)

#Pre-filtering dataset
#Removing rows with no counts.
#In this case we are only removing 0 counts as Michael Love (author) mentions in the [support page](#https://support.bioconductor.org/p/65256/) of Bioconductor that prefiltering of lowly expressed genes is usually not necessary when using DESeq2 as the pipeline includes independent filtering steps later and the presence of lowly expressed genes in this point does not affects results at the end of the workflow.
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)

#Normalization (transformation of counts)
#DESeq2 offers two transformation methods, variance stabilizing transformation (VST) and regularized-logarithm transformation (rlog). It is important to mention that this step is not needed for the gene expression analysis as DESeq2 works better on raw counts.

#We will test both transformation, it is important to consider that VST is computationally more efficient in larger datasets like the one used.
#VST
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
#rlog
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

#Scatterplots
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "Unnormalized"),
  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "Normalized - VST"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "Normalized - rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("Unnormalized", "Normalized - VST", "Normalized - rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)

#Sample distances-assessing similarity between samples

#Here we will generate a PCA plot and an MDS plot
#The objective is to analyze similarity between samples to determine if they fit the experimental design. 

#Multidimensional scaling analysis (MDS plot)
plotPCA(vsd, intgroup = c("dex", "cell"))

#MDS plot
#Multidimensional scaling (MDS) function in base R.
sampleDists <- dist(t(assay(vsd)))
sampleDists2 <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
sampleDistMatrix2 <- as.matrix(sampleDists2)
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
mds2 <- as.data.frame(colData(rld))  %>%
  cbind(cmdscale(sampleDistMatrix2))

ggplot(mds, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS-VST transformed")
ggplot(mds2, aes(x = `1`, y = `2`, color = dex, shape = cell)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS-rlog transformed")

###Code Section 2 - Main Analysis and Visualizations 

#Differential expression analysis main function generating a DESeqDataSet object.
#Is important to mention that we already defined the experimental design in a previous step as design = ~ cell + dex.
#This function performs 3 main steps: 1. Estimation of size factors. 2. Estimation of dispersion values for each gene. 3. Fitting a generalized linear model.
dds <- DESeq(dds)
class(dds)

#Generating results of the model indicating: base means across samples, log2 fold changes, standard errors, test statistics, p-values and adjusted p-values.
res <- results(dds)
#In the following plot, baseMean shows the average of the normalized count values. Log2FoldChange (lfc) shows how much the gene’s expression seems to have changed due to treatment with dexamethasone in comparison to untreated samples and lfcSE indicates the standard error of the lfc.
head(res)
summary(res)

#With the previous summary we can see that many genes have a response higher (2373) or lower (1949) than 0. Instead we want genes with a high response to DEX so a technique would be by raising log2 fold change threshold. The default was 0, by using 1 we are only including genes doubling or halving their response to DEX. This process is called lowering the false discovery rate.
resLFC1 <- results(dds, lfcThreshold=1)
summary(resLFC1)
#Now we have 138 genes (LFC >1) and 73 genes(LFC <1)

#Plotting results of gene expression

#Finding the gene with highest response
topGene <- rownames(resLFC1)[which.min(res$padj)]
#After a quick search in Ensembl database we find the gene is Monoamine oxidase A
#Plotting counts of Monoamine oxidase A gene as a response to DEX 
plotCounts(dds, gene = topGene, intgroup=c("dex"), main = "Counts of Monoamine oxidase A gene as a response to DEX")

#Bland–Altman (BA) plot
#Shows an overview for the distribution of the estimated coefficients in the model.
library("apeglm")
resultsNames(dds)
resLFC1 <- lfcShrink(dds, coef="dex_trt_vs_untrt", type="apeglm")
plotMA(resLFC1, ylim = c(-5, 5), main = "MA-plot of changes induced by treatment")
#Plot shows distribution of upregulated, downregulated (blue) and lowly expressed (grey) counts in data 
