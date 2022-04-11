#' ---
#' title: "FT CRISPR Biological QA"
#' author: "Nicolas Delhomme & Domenique Andre"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' ```{r no warn,echo=FALSE}
#' options(warn=-1)
#' ```

#' * Libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' * Helper functions
source(here("UPSCb-common/src/R/plot.multidensity.R"))
source(here("UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
pal <- brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Metadata
#' Sample information
samples <- read_delim(here("doc/Samples_correct2.txt"),
                      delim="\t",col_types=cols(
                        col_factor(),
                        col_factor(),
                        col_double(),
                        col_double()
                      ))

#' Extract some more info from the UserID
# samples <- samples %>% 
# samples %<>% 
samples %<>% mutate(Genotype=factor(sub("_.*","",UserID)),
                    Treatment=sub(".*_","",sub("_[0-9]$","",UserID)))

#' Some more manipulation to correct for T0/WT
samples$Treatment[samples$Treatment %in% samples$Genotype] <- "WT"
samples %<>% mutate(Genotype=relevel(Genotype,"T89")) %>% 
  mutate(Treatment=relevel(factor(Treatment),"WT"))


#' tx2gene translation table
tx2gene <- suppressMessages(read_delim(here("reference/annotation/tx2gene.tsv"),delim="\t",
                                       col_names=c("TXID","GENE")))

#' # Raw data
filelist <- list.files(here("data/salmonPotra02"), 
                       recursive = TRUE, 
                       pattern = "quant.sf",
                       full.names = TRUE)

#' Sanity check to ensure that the data is sorted according to the sample info
stopifnot(all(match(sub("_S.*","",basename(dirname(filelist))),samples$NGI_ID) == 1:nrow(samples)))

#' name the file list vector
names(filelist) <- samples$NGI_ID

#' Read the expression at the gene level
counts <- suppressMessages(round(tximport(files = filelist, 
                                          type = "salmon",
                                          tx2gene=tx2gene)$counts))

#' ## Quality Control
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' check how many genes are never expressed in a given condition
conditions=samples$Treatment
colSums(sapply(split.data.frame(t(counts!=0),conditions),colSums)==0)

save(counts, conditions, file = "counts_conditions.rda")

#' Let us take a look at the sequencing depth, colouring first by Genotype and 
#' then Treatment 
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

#' * Genotype
#' Some samples are more deeply sequenced than others but these are randomly
#' distributed across Genotype
ggplot(dat,aes(x,y,fill=Genotype)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' * Treatment
#' Neither by Treatment really
ggplot(dat,aes(x,y,fill=Treatment)) + geom_col() + 
  scale_y_continuous(name="reads") +
  theme(axis.text.x=element_text(angle=90,size=4),axis.title.x=element_blank())

#' Display the per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative gene coverage is as expected
ggplot(melt(log10(rowMeans(counts))),aes(x=value)) + 
  geom_density() + ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)")

#' The same is done for the individual samples colored by condition. The gene coverage 
#' across samples is extremely similar
dat <- as.data.frame(log10(counts)) %>% utils::stack() %>% 
  mutate(Genotype=samples$Genotype[match(ind,samples$NGI_ID)]) %>% 
  mutate(Treatment=samples$Treatment[match(ind,samples$NGI_ID)])

#' Color by Genotype
ggplot(dat,aes(x=values,group=ind,col=Genotype)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' Color by Treatment 
ggplot(dat,aes(x=values,group=ind,col=Treatment)) + 
  geom_density() + ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)")

#' ## Export
dir.create(here("data/analysis/salmon"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("data/analysis/salmon/raw-unormalised-gene-expression_data_corrected2.csv"))

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#' 
#' Note that the _design_ here is just expression as a function of `Genotype`. There is not
#' every combination of `Treatment` for every `Genotype`, so we cannot use a more complex model.
#' This is irrelevant here, as we are only doing the quality assessment of the data, but we 
#' need to keep it in mind for the differential expression analyses

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = samples,
  design = ~ Genotype)

#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
#' The sequencing depth is relatively variable (51% to 312 %), but that is due to
#' some samples being very deeply sequenced.
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' * Validation
#' 
#' The variance stabilisation worked adequately
#' 
meanSdPlot(vst[rowSums(counts)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)

#' * Genotype
#' The 2 first dimension explain a lot of the data variance
#' 
#' The Genotype 808 and 809 are from a separate experiment and cluster
#' together
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Genotype)],pch=19)
legend("topleft",
       fill=pal[1:nlevels(samples$Genotype)],
       legend=levels(samples$Genotype))

#' * Treatment
#' There are 3 groups (apart from the 808 and 809 Genotype) and
#' some of them have mixed Treatment. This is suspicious.
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Treatment)],pch=19)
legend("topleft",
       fill=pal[1:nlevels(samples$Treatment)],
       legend=levels(samples$Treatment))

par(mar=mar)

#' Check distribution of PCAs
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    samples)

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=Treatment,shape=Genotype,text=NGI_ID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' Is there a more general batch effect
pc.dat %<>% mutate(Batch=substr(samples$NGI_ID,8,8))
ggplot(pc.dat,aes(x=PC1,y=PC2,col=Batch)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts") +
  scale_x_continuous(paste("PC1 (",percent[1],"%)",sep="")) +
  scale_y_continuous(paste("PC2 (",percent[2],"%)",sep=""))

#' ### Heatmap
#' Filter for noise
#' A cutoff at a VST value of 5 leaves about 14000 genes, which is adequate for the QA
conds <- factor(paste(samples$Genotype,samples$Treatment))
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3)

#' * Heatmap of "all" genes
#' Taking into account all the genes (above a noise thresholds), the samples cluster
#' according to what we also see in the PCA, _i.e._ there is, depending on the time points, 
#' either a batch effect or some outliers 
hm <- heatmap.2(t(scale(t(vst[sels[[6]],]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,cexCol=.6,
                col=hpal)

#' Plot the sample dendrogram for better visualisation
plot(as.hclust(hm$colDendrogram),xlab="",sub="",main="Samples dendrogram",
     labels=paste(samples$Genotype,samples$Treatment),cex=0.6)

plot(as.hclust(hm$colDendrogram),xlab="",sub="",main="Samples dendrogram",
     labels=samples$UserID,cex=0.6)

plot(as.hclust(hm$colDendrogram),xlab="",sub="",main="Samples dendrogram",
     labels=substr(samples$NGI_ID,8,8),cex=0.6)


#' ## Conclusion
#' The data quality is overall good, the sequencing depth about constant.
#' 
#' Now remove outliers and save dds files for analysis.
#' Remove outliers and FT2 CRISPR data (by NGI_ID)
outliers <- c("P12108_101", "P12108_221", "P12108_231", #outliers
              "P12108_107", "P12108_108", "P12108_109", "P12108_116", "P12108_117", "P12108_118", #736 B1 LD, SDW15
              "P12108_125", "P12108_126", "P12108_127", "P12108_134", "P12108_135", "P12108_136", #736 B1 CTW2, W4
              "P12108_143", "P12108_144", "P12108_145", "P12108_152", "P12108_153", "P12108_154", #736 B1 W8, LDD7
              "P12108_207", "P12108_208", "P12108_209", "P12108_216", "P12108_217", "P12108_218", #736 B2 LD, SDW15
              "P12108_225", "P12108_226", "P12108_227", "P12108_234", "P12108_235", "P12108_236", #736 B2 CTW2, W4
              "P12108_243", "P12108_244", "P12108_245", "P12108_252", "P12108_253", "P12108_254", #736 B2 W8, LDD7
              "P12108_155", "P12108_156", "P12108_157", "P12108_158", "P12108_159", "P12108_160", "P12108_161", "P12108_162", "P12108_162", #FT2 B1
              "P12108_255", "P12108_256", "P12108_257", "P12108_258", "P12108_259", "P12108_260", "P12108_261", "P12108_262", "P12108_262") #FT2 B2
dds <- dds[,!dds$NGI_ID %in% outliers]


save(dds,file=here("data/analysis/salmon/dds_corrected.rda"))



#' 
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#'
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
