#' ---
#' title: "DESeq Differential Expression T89 vs FT1 CRISPR"
#' author: "Nicolas Delhomme & Domenique Andre"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' * Working directory

#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(plotly)
  library(RColorBrewer)
  library(tidyverse)
  library(VennDiagram)
})

#' * Helper files
suppressMessages(source(here("UPSCb-common/src/R/gopher.R")))
suppressMessages(source(here("UPSCb-common/src/R/plotMA.R")))
suppressMessages(source(here("UPSCb-common/src/R/volcanoPlot.R")))
suppressMessages(source(here("UPSCb-common/src/R/featureSelection.R")))

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds,vst,gene_id){
  sel <- grepl(gene_id,rownames(vst))
  stopifnot(sum(sel)==1)
  
  return(
    ggplot(bind_cols(as.data.frame(colData(dds)),
                     data.frame(value=vst[sel,])),
           aes(x=Treatment,y=value,col=Genotype,group=Genotype,text=NGI_ID)) +
      geom_point() + geom_smooth() +
      scale_y_continuous(name="VST expression") + 
      ggtitle(label=paste("Expression for: ",gene_id))
  )
}

"box_plot_LD" <- function(dds,vst,gene_id,treatment="LD"){
  row.sel <- grepl(gene_id,rownames(vst))
  col.sel <- dds$Treatment == treatment
  stopifnot(sum(row.sel)==1)
  
  return(
    ggplot(bind_cols(as.data.frame(colData(dds[,col.sel])),
                     data.frame(value=vst[row.sel,col.sel])),
           aes(x=Genotype,y=value,col=Genotype,group=Genotype,text=NGI_ID)) +
      geom_point() + geom_boxplot() +
      scale_y_continuous(name="VST expression") + 
      ggtitle(label=paste("Expression for: ",gene_id))
  )
}

"box_plot_SDW15" <- function(dds,vst,gene_id,treatment="SDW15"){
  row.sel <- grepl(gene_id,rownames(vst))
  col.sel <- dds$Treatment == treatment
  stopifnot(sum(row.sel)==1)
  
  return(
    ggplot(bind_cols(as.data.frame(colData(dds[,col.sel])),
                     data.frame(value=vst[row.sel,col.sel])),
           aes(x=Genotype,y=value,col=Genotype,group=Genotype,text=NGI_ID)) +
      geom_point() + geom_boxplot() +
      scale_y_continuous(name="VST expression") + 
      ggtitle(label=paste("Expression for: ",gene_id))
  )
}

"box_plot_CTW8" <- function(dds,vst,gene_id,treatment="CTW8"){
  row.sel <- grepl(gene_id,rownames(vst))
  col.sel <- dds$Treatment == treatment
  stopifnot(sum(row.sel)==1)
  
  return(
    ggplot(bind_cols(as.data.frame(colData(dds[,col.sel])),
                     data.frame(value=vst[row.sel,col.sel])),
           aes(x=Genotype,y=value,col=Genotype,group=Genotype,text=NGI_ID)) +
      geom_point() + geom_boxplot() +
      scale_y_continuous(name="VST expression") + 
      ggtitle(label=paste("Expression for: ",gene_id))
  )
}

"box_plot_LDD7" <- function(dds,vst,gene_id,treatment="LDD7"){
  row.sel <- grepl(gene_id,rownames(vst))
  col.sel <- dds$Treatment == treatment
  stopifnot(sum(row.sel)==1)
  
  return(
    ggplot(bind_cols(as.data.frame(colData(dds[,col.sel])),
                     data.frame(value=vst[row.sel,col.sel])),
           aes(x=Genotype,y=value,col=Genotype,group=Genotype,text=NGI_ID)) +
      geom_point() + geom_boxplot() +
      scale_y_continuous(name="VST expression") + 
      ggtitle(label=paste("Expression for: ",gene_id))
  )
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds)){
  
  if(length(contrast)==1){
    res <- results(dds,name=contrast)
  } else {
    res <- results(dds,contrast=contrast)
  }
  
  if(plot){
    par(mar=c(5,5,5,5))
    volcanoPlot(res)
    par(mar=mar)
  }
  
  sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj)
  
  if(verbose){
    message(sprintf("There are %s genes that are DE",sum(sel)))
  }
  
  if(export){
    if(!dir.exists(default_dir)){
      dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
    }
    write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
    write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
  }
  if(plot){
    if (sum(sel)>1){
      heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                distfun = pearson.dist,
                hclustfun = function(X){hclust(X,method="ward.D2")},
                trace="none",col=hpal,labRow = FALSE,
                labCol=labels[sample_sel], cexCol = 0.8, margins = c(7.1, 0.1))
    }
  }
  return(list(all=rownames(res[sel,]),
              up=rownames(res[sel & res$log2FoldChange>0,]),
              dn=rownames(res[sel & res$log2FoldChange<0,])))
}

#' # Data
load(here("data/analysis/salmon/dds_corrected.rda"))

#' * Remove the WT data from FT2 CRISPR (Treatment label "WT")
dds <- dds[,!colData(dds)$Treatment == "WT"]

#' * And reset the categorical variables (remove the missing levels)
dds$Treatment <- factor(dds$Treatment,levels=c("LD", "SDW15", "CTW2","CTW4","CTW8", "LDD7"))
dds$UserID <- factor(dds$UserID,levels=c("T89_LD_1", "T89_LD_2", "T89_LD_3", "T89_LD_4", "T89_LD_5", "T89_LD_6",
                                         "FT1_LD_1", "FT1_LD_2", "FT1_LD_3", "FT1_LD_4", "FT1_LD_5", "FT1_LD_6",
                                         "T89_SDW15_1", "T89_SDW15_2", "T89_SDW15_3", "T89_SDW15_4", "T89_SDW15_5", "T89_SDW15_6",
                                         "FT1_SDW15_1", "FT1_SDW15_2", "FT1_SDW15_3", "FT1_SDW15_4", "FT1_SDW15_5", "FT1_SDW15_6",
                                         "T89_CTW2_1", "T89_CTW2_2", "T89_CTW2_3", "T89_CTW2_4", "T89_CTW2_5", "T89_CTW2_6",
                                         "FT1_CTW2_1", "FT1_CTW2_2", "FT1_CTW2_3", "FT1_CTW2_4", "FT1_CTW2_5", "FT1_CTW2_6",
                                         "T89_CTW4_1", "T89_CTW4_2", "T89_CTW4_3", "T89_CTW4_4", "T89_CTW4_5", "T89_CTW4_6",
                                         "FT1_CTW4_1", "FT1_CTW4_2", "FT1_CTW4_3", "FT1_CTW4_4", "FT1_CTW4_5", "FT1_CTW4_6",
                                         "T89_CTW8_1", "T89_CTW8_2", "T89_CTW8_3", "T89_CTW8_4", "T89_CTW8_5", "T89_CTW8_6",
                                         "FT1_CTW8_1", "FT1_CTW8_2", "FT1_CTW8_3", "FT1_CTW8_4", "FT1_CTW8_5", "FT1_CTW8_6",
                                         "T89_LDD7_1", "T89_LDD7_2", "T89_LDD7_3", "T89_LDD7_4", "T89_LDD7_5", "T89_LDD7_6",
                                         "FT1_LDD7_1", "FT1_LDD7_2", "FT1_LDD7_3", "FT1_LDD7_4", "FT1_LDD7_5", "FT1_LDD7_6",
                                         ))
#' 
dds$Genotype <- droplevels(dds$Genotype)
dds$Batch <- factor(substr(dds$NGI_ID,8,8))

#' Finally update the design
design(dds) <- ~Genotype * Treatment

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

# save(vst,file=here("analysis/FT1_CRISPR_V2_app/vst_corrected.rda"))
# save(dds,file=here("analysis/FT1_CRISPR_V2_app/dds_corrected.rda"))
load(here("analysis/FT1_CRISPR_V2_app/vst_corrected.rda"))
load(here("analysis/FT1_CRISPR_V2_app/dds_corrected.rda"))
# write.csv(vst,"analysis/FT1_CRISPR_V2_app/vst.csv", row.names = TRUE)

#' ## Gene of interest
#' remove line 736 and outliers
sel.T89.FT1 <- dds$Genotype != "736" & ! dds$NGI_ID %in% c("P12108_101", "P12108_221", "P12108_231")

#' Heatmap of "all" genes
#' Create new order (by time point)
sel2 <- order(dds$UserID)
# re-order the columns of the dds
dds <- dds[,sel2]
# and the vst
vst <- vst[,sel2]
# double checking
stopifnot(colnames(vst)==colnames(dds))

#' Same as in biological QA, but only with T89 and FT1 and without outliers
conds <- factor(paste(dds$Genotype,dds$Treatment))[sel.T89.FT1]
sels <- rangeFeatureSelect(counts=vst[,sel.T89.FT1],
                           conditions=conds,
                           nrep=3)

hm <- heatmap.2(t(scale(t(vst[sels[[6]],sel.T89.FT1]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,cexCol=.6,
                col=hpal)

#' PCA for all samples
pc.T89.FT1 <- prcomp(t(vst))
percent.T89.FT1 <- round(summary(pc.T89.FT1)$importance[2,]*100)
pc.dat.T89.FT1 <- bind_cols(PC1=pc.T89.FT1$x[,1],
                            PC2=pc.T89.FT1$x[,2],
                            as.data.frame(colData(dds)))
p.T89.FT1 <- ggplot(pc.dat.T89.FT1,aes(x=PC1,y=PC2,col=Treatment,shape=Genotype,text=UserID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")
ggplotly(p.T89.FT1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent.T89.FT1[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent.T89.FT1[2],"%)",sep="")))


#' line_plot for all time points
#' 
#' ### FT1 Potra2n8c17315
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n8c17315")

#' ### FT2a Potra2n10c20842
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n10c20842")

#' ### FT2b Potra2n10c20839
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n10c20839")

#' ### SVP Potra2n7c16552
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n7c16552")

#' ### LAP1 Potra2n10c21057
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n10c21057")

#' ### CO1 Potra2n17c30991
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n17c30991")

#' ### CO2 Potra2n4c9430
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n4c9430")

#' ### GI Potra2n5c11095
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n5c11095")

#' ### GIL Potra2n2c5842
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n2c5842")

#' ### CENL1 Potra2n4c10201
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n4c10201")

#' ### CENL2 Potra2n9c18678
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n9c18678")

#' ### GA20ox-1 Potra2n12c24804
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n12c24804")

#' ### GA20ox-2 Potra2n15c28070
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n15c28070")

#' ### GA2 oxidase Potra2n14c27286
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n14c27286")

#' ### PhyA Potra2n13c26308
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n13c26308")

#' ### PhyB1 Potra2n8c17574
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n8c17574")

#' ### PhyB2 Potra2n10c21137
line_plot(dds[,sel.T89.FT1],vst[,sel.T89.FT1],"Potra2n10c21137")

#' # Differential Expression
#' 
#' According to the biological QA, the only outliers are 101 and 221
#' 
#' ## T89 time series, compare two consecutive time points
#' ### * T89 - LD _vs._ SDW15
#' outliers: P12108_101
sel.T89.LD.SDW15 <- dds$Genotype=="T89" & dds$Treatment %in% c("LD","SDW15") & dds$NGI_ID != "P12108_101"
ddsEx.T89.LD.SDW15 <- dds[,sel.T89.LD.SDW15]
design(ddsEx.T89.LD.SDW15) <- ~Treatment
ddsEx.T89.LD.SDW15$Treatment <- droplevels(ddsEx.T89.LD.SDW15$Treatment)
ddsEx.T89.LD.SDW15 <- DESeq(ddsEx.T89.LD.SDW15)
plotDispEsts(ddsEx.T89.LD.SDW15)
resultsNames(ddsEx.T89.LD.SDW15)
T89.LD.SDW15 <- extract_results(ddsEx.T89.LD.SDW15,vst[,sel.T89.LD.SDW15],"Treatment_SDW15_vs_LD",
                                default_prefix="T89_SDW15-vs-LD_",
                                labels=ddsEx.T89.LD.SDW15$UserID)
barplot(sapply(T89.LD.SDW15,length))


#' ### * T89 - SDW15 _vs._ CTW2
#' outliers: P12108_221
sel.T89.SDW15.CTW2 <- dds$Genotype=="T89" & dds$Treatment %in% c("SDW15","CTW2") & dds$NGI_ID != "P12108_221"
ddsEx.T89.SDW15.CTW2 <- dds[,sel.T89.SDW15.CTW2]
design(ddsEx.T89.SDW15.CTW2) <- ~Treatment
ddsEx.T89.SDW15.CTW2$Treatment <- droplevels(ddsEx.T89.SDW15.CTW2$Treatment)
ddsEx.T89.SDW15.CTW2 <- DESeq(ddsEx.T89.SDW15.CTW2)
plotDispEsts(ddsEx.T89.SDW15.CTW2)
resultsNames(ddsEx.T89.SDW15.CTW2)
T89.SDW15.CTW2 <- extract_results(ddsEx.T89.SDW15.CTW2,vst[,sel.T89.SDW15.CTW2],"Treatment_CTW2_vs_SDW15",
                                  default_prefix="T89_CTW2-vs-SDW15_",
                                  labels=ddsEx.T89.SDW15.CTW2$UserID)
barplot(sapply(T89.SDW15.CTW2,length))

#' ### * T89 - CTW2 _vs._ CTW4
#' outliers: P12108_221
sel.T89.CTW2.CTW4 <- dds$Genotype=="T89" & dds$Treatment %in% c("CTW2","CTW4") & dds$NGI_ID != "P12108_221"
ddsEx.T89.CTW2.CTW4 <- dds[,sel.T89.CTW2.CTW4]
design(ddsEx.T89.CTW2.CTW4) <- ~Treatment
ddsEx.T89.CTW2.CTW4$Treatment <- droplevels(ddsEx.T89.CTW2.CTW4$Treatment)
ddsEx.T89.CTW2.CTW4 <- DESeq(ddsEx.T89.CTW2.CTW4)
plotDispEsts(ddsEx.T89.CTW2.CTW4)
resultsNames(ddsEx.T89.CTW2.CTW4)
T89.CTW2.CTW4 <- extract_results(ddsEx.T89.CTW2.CTW4,vst[,sel.T89.CTW2.CTW4],"Treatment_CTW4_vs_CTW2",
                                 default_prefix="T89_CTW4-vs-CTW2_",
                                 labels=ddsEx.T89.CTW2.CTW4$UserID)
barplot(sapply(T89.CTW2.CTW4,length))

#' ### * T89 - CTW4 _vs._ CTW8
#' outliers: 
sel.T89.CTW4.CTW8 <- dds$Genotype=="T89" & dds$Treatment %in% c("CTW4","CTW8")
ddsEx.T89.CTW4.CTW8 <- dds[,sel.T89.CTW4.CTW8]
design(ddsEx.T89.CTW4.CTW8) <- ~Treatment
ddsEx.T89.CTW4.CTW8$Treatment <- droplevels(ddsEx.T89.CTW4.CTW8$Treatment)
ddsEx.T89.CTW4.CTW8 <- DESeq(ddsEx.T89.CTW4.CTW8)
plotDispEsts(ddsEx.T89.CTW4.CTW8)
resultsNames(ddsEx.T89.CTW4.CTW8)
T89.CTW4.CTW8 <- extract_results(ddsEx.T89.CTW4.CTW8,vst[,sel.T89.CTW4.CTW8],"Treatment_CTW8_vs_CTW4",
                                 default_prefix="T89_CTW8-vs-CTW4_",
                                 labels=ddsEx.T89.CTW4.CTW8$UserID)
barplot(sapply(T89.CTW4.CTW8,length))

#' ### * T89 - CTW8 _vs._ LDD7
#' outliers:
sel.T89.CTW8.LDD7 <- dds$Genotype=="T89" & dds$Treatment %in% c("CTW8","LDD7")
ddsEx.T89.CTW8.LDD7 <- dds[,sel.T89.CTW8.LDD7]
design(ddsEx.T89.CTW8.LDD7) <- ~Treatment
ddsEx.T89.CTW8.LDD7$Treatment <- droplevels(ddsEx.T89.CTW8.LDD7$Treatment)
ddsEx.T89.CTW8.LDD7 <- DESeq(ddsEx.T89.CTW8.LDD7)
plotDispEsts(ddsEx.T89.CTW8.LDD7)
resultsNames(ddsEx.T89.CTW8.LDD7)
T89.CTW8.LDD7 <- extract_results(ddsEx.T89.CTW8.LDD7,vst[,sel.T89.CTW8.LDD7],"Treatment_LDD7_vs_CTW8",
                                 default_prefix="T89_LDD7-vs-CTW8_",
                                 labels=ddsEx.T89.CTW8.LDD7$UserID)
barplot(sapply(T89.CTW8.LDD7,length))

#' ## Gene expression over time
#' Barplot to show number of DE genes from timepoint to timepoint in T89
barplot(sapply(list(LD_SDW15=T89.LD.SDW15,SDW15_CTW2=T89.SDW15.CTW2,
                    CTW2_CTW4=T89.CTW2.CTW4,CTW4_CTW8=T89.CTW4.CTW8,
                    CTW8_LDD7=T89.CTW8.LDD7),
               elementNROWS),beside=TRUE,legend.text = TRUE)

#' Venn diagram for all DE genes at different time points
#' all DE genes
grid.newpage()
grid.draw(venn.diagram(list(LD_SDW15=T89.LD.SDW15$all, SDW15_CTW2=T89.SDW15.CTW2$all, 
                            CTW2_CTW4=T89.CTW2.CTW4$all, CTW4_CTW8=T89.CTW4.CTW8$all, 
                            CTW8_LDD7=T89.CTW8.LDD7$all),filename = NULL,fill=pal[1:5]))

grid.newpage()
grid.draw(venn.diagram(list(SDW15_CTW2=T89.SDW15.CTW2$all, 
                            CTW2_CTW4=T89.CTW2.CTW4$all, CTW4_CTW8=T89.CTW4.CTW8$all, 
                            CTW8_LDD7=T89.CTW8.LDD7$all),filename = NULL,fill=pal[1:4]))
#' up regulated genes
grid.newpage()
grid.draw(venn.diagram(list(LD_SDW15=T89.LD.SDW15$up, SDW15_CTW2=T89.SDW15.CTW2$up, 
                            CTW2_CTW4=T89.CTW2.CTW4$up, CTW4_CTW8=T89.CTW4.CTW8$up, 
                            CTW8_LDD7=T89.CTW8.LDD7$up),filename = NULL,fill=pal[1:5]))

grid.newpage()
grid.draw(venn.diagram(list(SDW15_CTW2=T89.SDW15.CTW2$up, 
                            CTW2_CTW4=T89.CTW2.CTW4$up, CTW4_CTW8=T89.CTW4.CTW8$up, 
                            CTW8_LDD7=T89.CTW8.LDD7$up),filename = NULL,fill=pal[1:4]))
#' down regulated genes
grid.newpage()
grid.draw(venn.diagram(list(LD_SDW15=T89.LD.SDW15$dn, SDW15_CTW2=T89.SDW15.CTW2$dn, 
                            CTW2_CTW4=T89.CTW2.CTW4$dn, CTW4_CTW8=T89.CTW4.CTW8$dn, 
                            CTW8_LDD7=T89.CTW8.LDD7$dn),filename = NULL,fill=pal[1:5]))


#' shared genes between T89 CTW2 to CTW4 to CTW8
grid.newpage()
grid.draw(venn.diagram(list(CTW2_CTW4=T89.CTW2.CTW4$all, CTW4_CTW8=T89.CTW4.CTW8$all),filename = NULL,fill=pal[1:2]))
grid.newpage()
grid.draw(venn.diagram(list(CTW2_CTW4=T89.CTW2.CTW4$up, CTW4_CTW8=T89.CTW4.CTW8$dn),filename = NULL,fill=pal[1:2]))
grid.newpage()
grid.draw(venn.diagram(list(CTW2_CTW4=T89.CTW2.CTW4$dn, CTW4_CTW8=T89.CTW4.CTW8$up),filename = NULL,fill=pal[1:2]))


#' ### * T89 - SDW15 _vs._ CTW8
#' check the difference between ecodormancy before and after endodormancy
sel.T89.SDW15.CTW8 <- dds$Genotype=="T89" & dds$Treatment %in% c("SDW15","CTW8")
ddsEx.T89.SDW15.CTW8 <- dds[,sel.T89.SDW15.CTW8]
design(ddsEx.T89.SDW15.CTW8) <- ~Treatment
ddsEx.T89.SDW15.CTW8$Treatment <- droplevels(ddsEx.T89.SDW15.CTW8$Treatment)
ddsEx.T89.SDW15.CTW8 <- DESeq(ddsEx.T89.SDW15.CTW8)
plotDispEsts(ddsEx.T89.SDW15.CTW8)
resultsNames(ddsEx.T89.SDW15.CTW8)
T89.SDW15.CTW8 <- extract_results(ddsEx.T89.SDW15.CTW8,vst[,sel.T89.SDW15.CTW8],"Treatment_CTW8_vs_SDW15",
                                  default_prefix="T89_CTW8-vs-SDW15_",
                                  labels=ddsEx.T89.SDW15.CTW8$UserID)
barplot(sapply(T89.SDW15.CTW8,length))

sel.T89.CTW2.CTW8 <- dds$Genotype=="T89" & dds$Treatment %in% c("CTW2","CTW8")
ddsEx.T89.CTW2.CTW8 <- dds[,sel.T89.CTW2.CTW8]
design(ddsEx.T89.CTW2.CTW8) <- ~Treatment
ddsEx.T89.CTW2.CTW8$Treatment <- droplevels(ddsEx.T89.CTW2.CTW8$Treatment)
ddsEx.T89.CTW2.CTW8 <- DESeq(ddsEx.T89.CTW2.CTW8)
plotDispEsts(ddsEx.T89.CTW2.CTW8)
resultsNames(ddsEx.T89.CTW2.CTW8)
T89.CTW2.CTW8 <- extract_results(ddsEx.T89.CTW2.CTW8,vst[,sel.T89.CTW2.CTW8],"Treatment_CTW8_vs_CTW2",
                                 default_prefix="T89_CTW8-vs-CTW2_",
                                 labels=ddsEx.T89.CTW2.CTW8$UserID)
barplot(sapply(T89.CTW2.CTW8,length))

sel.T89.LD.CTW8 <- dds$Genotype=="T89" & dds$Treatment %in% c("LD","CTW8")
ddsEx.T89.LD.CTW8 <- dds[,sel.T89.LD.CTW8]
design(ddsEx.T89.LD.CTW8) <- ~Treatment
ddsEx.T89.LD.CTW8$Treatment <- droplevels(ddsEx.T89.LD.CTW8$Treatment)
ddsEx.T89.LD.CTW8 <- DESeq(ddsEx.T89.LD.CTW8)
plotDispEsts(ddsEx.T89.LD.CTW8)
resultsNames(ddsEx.T89.LD.CTW8)
T89.LD.CTW8 <- extract_results(ddsEx.T89.LD.CTW8,vst[,sel.T89.LD.CTW8],"Treatment_CTW8_vs_LD",
                               default_prefix="T89_CTW8-vs-LD_",
                               labels=ddsEx.T89.LD.CTW8$UserID)
barplot(sapply(T89.LD.CTW8,length))


#' ## FT1 time series, compare two consecutive time points
#' ### * FT1 - LD _vs._ SDW15
#' outliers:
sel.FT1.LD.SDW15 <- dds$Genotype=="FT1" & dds$Treatment %in% c("LD","SDW15")
ddsEx.FT1.LD.SDW15 <- dds[,sel.FT1.LD.SDW15]
design(ddsEx.FT1.LD.SDW15) <- ~Treatment
ddsEx.FT1.LD.SDW15$Treatment <- droplevels(ddsEx.FT1.LD.SDW15$Treatment)
ddsEx.FT1.LD.SDW15 <- DESeq(ddsEx.FT1.LD.SDW15)
plotDispEsts(ddsEx.FT1.LD.SDW15)
resultsNames(ddsEx.FT1.LD.SDW15)
FT1.LD.SDW15 <- extract_results(ddsEx.FT1.LD.SDW15,vst[,sel.FT1.LD.SDW15],"Treatment_SDW15_vs_LD",
                                default_prefix="FT1_SDW15-vs-LD_",
                                labels=ddsEx.FT1.LD.SDW15$UserID)
barplot(sapply(FT1.LD.SDW15,length))


#' ### * FT1 - SDW15 _vs._ CTW2
#' outliers: 
sel.FT1.SDW15.CTW2 <- dds$Genotype=="FT1" & dds$Treatment %in% c("SDW15","CTW2")
ddsEx.FT1.SDW15.CTW2 <- dds[,sel.FT1.SDW15.CTW2]
design(ddsEx.FT1.SDW15.CTW2) <- ~Treatment
ddsEx.FT1.SDW15.CTW2$Treatment <- droplevels(ddsEx.FT1.SDW15.CTW2$Treatment)
ddsEx.FT1.SDW15.CTW2 <- DESeq(ddsEx.FT1.SDW15.CTW2)
plotDispEsts(ddsEx.FT1.SDW15.CTW2)
resultsNames(ddsEx.FT1.SDW15.CTW2)
FT1.SDW15.CTW2 <- extract_results(ddsEx.FT1.SDW15.CTW2,vst[,sel.FT1.SDW15.CTW2],"Treatment_CTW2_vs_SDW15",
                                  default_prefix="FT1_CTW2-vs-SDW15_",
                                  labels=ddsEx.FT1.SDW15.CTW2$UserID)
barplot(sapply(FT1.SDW15.CTW2,length))

#' ### * FT1 - CTW2 _vs._ CTW4
#' outliers: 
sel.FT1.CTW2.CTW4 <- dds$Genotype=="FT1" & dds$Treatment %in% c("CTW2","CTW4")
ddsEx.FT1.CTW2.CTW4 <- dds[,sel.FT1.CTW2.CTW4]
design(ddsEx.FT1.CTW2.CTW4) <- ~Treatment
ddsEx.FT1.CTW2.CTW4$Treatment <- droplevels(ddsEx.FT1.CTW2.CTW4$Treatment)
ddsEx.FT1.CTW2.CTW4 <- DESeq(ddsEx.FT1.CTW2.CTW4)
plotDispEsts(ddsEx.FT1.CTW2.CTW4)
resultsNames(ddsEx.FT1.CTW2.CTW4)
FT1.CTW2.CTW4 <- extract_results(ddsEx.FT1.CTW2.CTW4,vst[,sel.FT1.CTW2.CTW4],"Treatment_CTW4_vs_CTW2",
                                 default_prefix="FT1_CTW4-vs-CTW2_",
                                 labels=ddsEx.FT1.CTW2.CTW4$UserID)
barplot(sapply(FT1.CTW2.CTW4,length))

#' ### * FT1 - CTW4 _vs._ CTW8
#' outliers: 
sel.FT1.CTW4.CTW8 <- dds$Genotype=="FT1" & dds$Treatment %in% c("CTW4","CTW8")
ddsEx.FT1.CTW4.CTW8 <- dds[,sel.FT1.CTW4.CTW8]
design(ddsEx.FT1.CTW4.CTW8) <- ~Treatment
ddsEx.FT1.CTW4.CTW8$Treatment <- droplevels(ddsEx.FT1.CTW4.CTW8$Treatment)
ddsEx.FT1.CTW4.CTW8 <- DESeq(ddsEx.FT1.CTW4.CTW8)
plotDispEsts(ddsEx.FT1.CTW4.CTW8)
resultsNames(ddsEx.FT1.CTW4.CTW8)
FT1.CTW4.CTW8 <- extract_results(ddsEx.FT1.CTW4.CTW8,vst[,sel.FT1.CTW4.CTW8],"Treatment_CTW8_vs_CTW4",
                                 default_prefix="FT1_CTW8-vs-CTW4_",
                                 labels=ddsEx.FT1.CTW4.CTW8$UserID)
barplot(sapply(FT1.CTW4.CTW8,length))

#' ### * FT1 - CTW8 _vs._ LDD7
#' outliers:
sel.FT1.CTW8.LDD7 <- dds$Genotype=="FT1" & dds$Treatment %in% c("CTW8","LDD7")
ddsEx.FT1.CTW8.LDD7 <- dds[,sel.FT1.CTW8.LDD7]
design(ddsEx.FT1.CTW8.LDD7) <- ~Treatment
ddsEx.FT1.CTW8.LDD7$Treatment <- droplevels(ddsEx.FT1.CTW8.LDD7$Treatment)
ddsEx.FT1.CTW8.LDD7 <- DESeq(ddsEx.FT1.CTW8.LDD7)
plotDispEsts(ddsEx.FT1.CTW8.LDD7)
resultsNames(ddsEx.FT1.CTW8.LDD7)
FT1.CTW8.LDD7 <- extract_results(ddsEx.FT1.CTW8.LDD7,vst[,sel.FT1.CTW8.LDD7],"Treatment_LDD7_vs_CTW8",
                                 default_prefix="FT1_LDD7-vs-CTW8_",
                                 labels=ddsEx.FT1.CTW8.LDD7$UserID)
barplot(sapply(FT1.CTW8.LDD7,length))

#' ## Gene expression over time
#' Barplot to show number of DE genes from timepoint to timepoint in FT1
barplot(sapply(list(LD_SDW15=FT1.LD.SDW15,SDW15_CTW2=FT1.SDW15.CTW2,
                    CTW2_CTW4=FT1.CTW2.CTW4,CTW4_CTW8=FT1.CTW4.CTW8,
                    CTW8_LDD7=FT1.CTW8.LDD7),
               elementNROWS),beside=TRUE,legend.text = TRUE)

#' Venn diagram for all DE genes at different time points
#' all DE genes
grid.newpage()
grid.draw(venn.diagram(list(LD_SDW15=FT1.LD.SDW15$all, SDW15_CTW2=FT1.SDW15.CTW2$all, CTW2_CTW4=FT1.CTW2.CTW4$all, CTW4_CTW8=FT1.CTW4.CTW8$all, CTW8_LDD7=FT1.CTW8.LDD7$all),filename = NULL,fill=pal[1:5]))
#' up regulated genes
grid.newpage()
grid.draw(venn.diagram(list(LD_SDW15=FT1.LD.SDW15$up, SDW15_CTW2=FT1.SDW15.CTW2$up, CTW2_CTW4=FT1.CTW2.CTW4$up, CTW4_CTW8=FT1.CTW4.CTW8$up, CTW8_LDD7=FT1.CTW8.LDD7$up),filename = NULL,fill=pal[1:5]))
#' down regulated genes
grid.newpage()
grid.draw(venn.diagram(list(LD_SDW15=FT1.LD.SDW15$dn, SDW15_CTW2=FT1.SDW15.CTW2$dn, CTW2_CTW4=FT1.CTW2.CTW4$dn, CTW4_CTW8=FT1.CTW4.CTW8$dn, CTW8_LDD7=FT1.CTW8.LDD7$dn),filename = NULL,fill=pal[1:5]))

#' shared genes between T89 and FT1 that change from CTW4 to CTW8
grid.newpage()
grid.draw(venn.diagram(list(T89=T89.CTW4.CTW8$all, FT1=FT1.CTW4.CTW8$all),filename = NULL,fill=pal[1:2]))


#' ## Compare T89 and FT1 at each time point. See when there start being differences.
#' ### * LD  T89 _vs._ FT1
#' outlier: P12108_101
sel.LD.T89.FT1 <- dds$Treatment=="LD" & dds$Genotype %in% c("T89","FT1") & dds$NGI_ID != "P12108_101"
ddsEx.LD.T89.FT1 <- dds[,sel.LD.T89.FT1]
design(ddsEx.LD.T89.FT1) <- ~Genotype
ddsEx.LD.T89.FT1$Genotype <- droplevels(ddsEx.LD.T89.FT1$Genotype)
ddsEx.LD.T89.FT1 <- DESeq(ddsEx.LD.T89.FT1)
plotDispEsts(ddsEx.LD.T89.FT1)
resultsNames(ddsEx.LD.T89.FT1)
LD.T89.FT1 <- extract_results(ddsEx.LD.T89.FT1,vst[,sel.LD.T89.FT1],"Genotype_FT1_vs_T89",
                              default_prefix="LD_FT1-vs-T89_",
                              labels=ddsEx.LD.T89.FT1$UserID)
barplot(sapply(LD.T89.FT1,length))

#' ### * SDW15  T89 _vs._ FT1
#' outlier: 
sel.SDW15.T89.FT1 <- dds$Treatment=="SDW15" & dds$Genotype %in% c("T89","FT1")
ddsEx.SDW15.T89.FT1 <- dds[,sel.SDW15.T89.FT1]
design(ddsEx.SDW15.T89.FT1) <- ~Genotype
ddsEx.SDW15.T89.FT1$Genotype <- droplevels(ddsEx.SDW15.T89.FT1$Genotype)
ddsEx.SDW15.T89.FT1 <- DESeq(ddsEx.SDW15.T89.FT1)
plotDispEsts(ddsEx.SDW15.T89.FT1)
resultsNames(ddsEx.SDW15.T89.FT1)
SDW15.T89.FT1 <- extract_results(ddsEx.SDW15.T89.FT1,vst[,sel.SDW15.T89.FT1],"Genotype_FT1_vs_T89",
                                 default_prefix="SDW15_FT1-vs-T89_",
                                 labels=ddsEx.SDW15.T89.FT1$UserID)
barplot(sapply(SDW15.T89.FT1,length))

#' ### * CTW2  T89 _vs._ FT1
#' outlier: P12108_221 
sel.CTW2.T89.FT1 <- dds$Treatment=="CTW2" & dds$Genotype %in% c("T89","FT1") & dds$NGI_ID != "P12108_221" 
ddsEx.CTW2.T89.FT1 <- dds[,sel.CTW2.T89.FT1]
design(ddsEx.CTW2.T89.FT1) <- ~Genotype
ddsEx.CTW2.T89.FT1$Genotype <- droplevels(ddsEx.CTW2.T89.FT1$Genotype)
ddsEx.CTW2.T89.FT1 <- DESeq(ddsEx.CTW2.T89.FT1)
plotDispEsts(ddsEx.CTW2.T89.FT1)
resultsNames(ddsEx.CTW2.T89.FT1)
CTW2.T89.FT1 <- extract_results(ddsEx.CTW2.T89.FT1,vst[,sel.CTW2.T89.FT1],"Genotype_FT1_vs_T89",
                                default_prefix="CTW2_FT1-vs-T89_",
                                labels=ddsEx.CTW2.T89.FT1$UserID)
barplot(sapply(CTW2.T89.FT1,length))

#' ### * CTW4  T89 _vs._ FT1
#' outlier: 
sel.CTW4.T89.FT1 <- dds$Treatment=="CTW4" & dds$Genotype %in% c("T89","FT1")
ddsEx.CTW4.T89.FT1 <- dds[,sel.CTW4.T89.FT1]
design(ddsEx.CTW4.T89.FT1) <- ~Genotype
ddsEx.CTW4.T89.FT1$Genotype <- droplevels(ddsEx.CTW4.T89.FT1$Genotype)
ddsEx.CTW4.T89.FT1 <- DESeq(ddsEx.CTW4.T89.FT1)
plotDispEsts(ddsEx.CTW4.T89.FT1)
resultsNames(ddsEx.CTW4.T89.FT1)
CTW4.T89.FT1 <- extract_results(ddsEx.CTW4.T89.FT1,vst[,sel.CTW4.T89.FT1],"Genotype_FT1_vs_T89",
                                default_prefix="CTW4_FT1-vs-T89_",
                                labels=ddsEx.CTW4.T89.FT1$UserID)
barplot(sapply(CTW4.T89.FT1,length))

#' ### * CTW8  T89 _vs._ FT1
#' outlier:
sel.CTW8.T89.FT1 <- dds$Treatment=="CTW8" & dds$Genotype %in% c("T89","FT1")
ddsEx.CTW8.T89.FT1 <- dds[,sel.CTW8.T89.FT1]
design(ddsEx.CTW8.T89.FT1) <- ~Genotype
ddsEx.CTW8.T89.FT1$Genotype <- droplevels(ddsEx.CTW8.T89.FT1$Genotype)
ddsEx.CTW8.T89.FT1 <- DESeq(ddsEx.CTW8.T89.FT1)
plotDispEsts(ddsEx.CTW8.T89.FT1)
resultsNames(ddsEx.CTW8.T89.FT1)
CTW8.T89.FT1 <- extract_results(ddsEx.CTW8.T89.FT1,vst[,sel.CTW8.T89.FT1],"Genotype_FT1_vs_T89",
                                default_prefix="CTW8_FT1-vs-T89_",
                                labels=ddsEx.CTW8.T89.FT1$UserID)
barplot(sapply(CTW8.T89.FT1,length))

#' ### * LDD7  T89 _vs._ FT1
#' outlier: 
sel.LDD7.T89.FT1 <- dds$Treatment=="LDD7" & dds$Genotype %in% c("T89","FT1")
ddsEx.LDD7.T89.FT1 <- dds[,sel.LDD7.T89.FT1]
design(ddsEx.LDD7.T89.FT1) <- ~Genotype
ddsEx.LDD7.T89.FT1$Genotype <- droplevels(ddsEx.LDD7.T89.FT1$Genotype)
ddsEx.LDD7.T89.FT1 <- DESeq(ddsEx.LDD7.T89.FT1)
plotDispEsts(ddsEx.LDD7.T89.FT1)
resultsNames(ddsEx.LDD7.T89.FT1)
LDD7.T89.FT1 <- extract_results(ddsEx.LDD7.T89.FT1,vst[,sel.LDD7.T89.FT1],"Genotype_FT1_vs_T89",
                                default_prefix="LDD7_FT1-vs-T89_",
                                labels=ddsEx.LDD7.T89.FT1$UserID)
barplot(sapply(LDD7.T89.FT1,length))

#' ## Gene expression over time for each time point between T89 and FT1 CRISPR
barplot(sapply(list(LD=LD.T89.FT1,SDW15=SDW15.T89.FT1,
                    CTW2=CTW2.T89.FT1,CTW4=CTW4.T89.FT1,
                    CTW8=CTW8.T89.FT1,LDD7=LDD7.T89.FT1),elementNROWS),
        beside=TRUE,legend.text = TRUE)

#' Venn diagram for all DE genes at different time points
#' all DE genes
grid.newpage()
grid.draw(venn.diagram(list(SDW15=SDW15.T89.FT1$all, CTW2=CTW2.T89.FT1$all, CTW4=CTW4.T89.FT1$all, 
                            CTW8=CTW8.T89.FT1$all, LDD7=LDD7.T89.FT1$all),filename = NULL,fill=pal[1:5]))
#' up regulated genes
grid.newpage()
grid.draw(venn.diagram(list(SDW15=SDW15.T89.FT1$up, CTW2=CTW2.T89.FT1$up, CTW4=CTW4.T89.FT1$up, 
                            CTW8=CTW8.T89.FT1$up, LDD7=LDD7.T89.FT1$up),filename = NULL,fill=pal[1:5]))
#' down regulated genes
grid.newpage()
grid.draw(venn.diagram(list(SDW15=SDW15.T89.FT1$dn, CTW2=CTW2.T89.FT1$dn, CTW4=CTW4.T89.FT1$dn, 
                            CTW8=CTW8.T89.FT1$dn, LDD7=LDD7.T89.FT1$dn),filename = NULL,fill=pal[1:5]))


#' Venn diagram for DE genes at T89 CTW4-CTW8 vs FT1 CTW4-CTW8
grid.newpage()
grid.draw(venn.diagram(list(T89=T89.CTW4.CTW8$all, FT1=FT1.CTW4.CTW8$all),filename = NULL,fill=pal[1:2]))
grid.newpage()
grid.draw(venn.diagram(list(T89=T89.CTW4.CTW8$up, FT1=FT1.CTW4.CTW8$up),filename = NULL,fill=pal[1:2]))
grid.newpage()
grid.draw(venn.diagram(list(T89=T89.CTW4.CTW8$dn, FT1=FT1.CTW4.CTW8$dn),filename = NULL,fill=pal[1:2]))
grid.newpage()
grid.draw(venn.diagram(list(T89=T89.CTW4.CTW8$up, FT1=FT1.CTW4.CTW8$dn),filename = NULL,fill=pal[1:2]))
grid.newpage()
grid.draw(venn.diagram(list(T89=T89.CTW4.CTW8$dn, FT1=FT1.CTW4.CTW8$up),filename = NULL,fill=pal[1:2]))


#' #' # GO Enrichment
#' #'
#' #' Here we perform a GO enrichment of the DE genes
#' #'
#' #' We use as background, the genes that have any level of expression in
#' #' the conditions we looked at. We could be more stringent and set some
#' #' noise filter (similarly to what we did in the Biological QA)
#' #'
#' #'Create a directory to store the GO output
#' dir.create(here("data/analysis/GO"),showWarnings = FALSE)
#' #' #'
#' #' ## GO enrichment for different time points in T89
#' 
#' ### GO enrichment T89 LD _vs._ SDW15
#' #' all DE genes
#' enr.T89.LD.SDW15 <- gopher(genes=T89.LD.SDW15$all, background=rownames(vst)[rowSums(counts(ddsEx.T89.LD.SDW15)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.LD.SDW15$go,file = here("data/analysis/GO/T89.LD.SDW15-GO-enrichment.csv"))
#' write_delim(enr.T89.LD.SDW15$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.LD.SDW15-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.T89.LD.SDW15.up <- gopher(genes=T89.LD.SDW15$up, background=rownames(vst)[rowSums(counts(ddsEx.T89.LD.SDW15)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.LD.SDW15.up$go,file = here("data/analysis/GO/T89.LD.SDW15.up-GO-enrichment.csv"))
#' write_delim(enr.T89.LD.SDW15.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.LD.SDW15.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.T89.LD.SDW15.dn <- gopher(genes=T89.LD.SDW15$dn, background=rownames(vst)[rowSums(counts(ddsEx.T89.LD.SDW15)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.LD.SDW15.dn$go,file = here("data/analysis/GO/T89.LD.SDW15-dn-GO-enrichment.csv"))
#' write_delim(enr.T89.LD.SDW15.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.LD.SDW15.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment T89 SDW15 _vs._ CTW2
#' #' all DE genes
#' enr.T89.SDW15.CTW2 <- gopher(genes=T89.SDW15.CTW2$all, background=rownames(vst)[rowSums(counts(ddsEx.T89.SDW15.CTW2)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.SDW15.CTW2$go,file = here("data/analysis/GO/T89.SDW15.CTW2-GO-enrichment.csv"))
#' write_delim(enr.T89.SDW15.CTW2$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.SDW15.CTW2-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.T89.SDW15.CTW2.up <- gopher(genes=T89.SDW15.CTW2$up, background=rownames(vst)[rowSums(counts(ddsEx.T89.SDW15.CTW2)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.SDW15.CTW2.up$go,file = here("data/analysis/GO/T89.SDW15.CTW2.up-GO-enrichment.csv"))
#' write_delim(enr.T89.SDW15.CTW2.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.SDW15.CTW2.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.T89.SDW15.CTW2.dn <- gopher(genes=T89.SDW15.CTW2$dn, background=rownames(vst)[rowSums(counts(ddsEx.T89.SDW15.CTW2)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.SDW15.CTW2.dn$go,file = here("data/analysis/GO/T89.SDW15.CTW2.dn-GO-enrichment.csv"))
#' write_delim(enr.T89.SDW15.CTW2.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.SDW15.CTW2.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment T89 CTW2 _vs._ CTW4
#' #' all DE genes
#' enr.T89.CTW2.CTW4 <- gopher(genes=T89.CTW2.CTW4$all, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW2.CTW4)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW2.CTW4$go,file = here("data/analysis/GO/T89.CTW2.CTW4-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW2.CTW4$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW2.CTW4-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.T89.CTW2.CTW4.up <- gopher(genes=T89.CTW2.CTW4$up, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW2.CTW4)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW2.CTW4.up$go,file = here("data/analysis/GO/T89.CTW2.CTW4.up-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW2.CTW4.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW2.CTW4.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.T89.CTW2.CTW4.dn <- gopher(genes=T89.CTW2.CTW4$dn, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW2.CTW4)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW2.CTW4.dn$go,file = here("data/analysis/GO/T89.CTW2.CTW4.dn-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW2.CTW4.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW2.CTW4.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment T89 CTW4 _vs._ CTW8
#' #' all DE genes
#' enr.T89.CTW4.CTW8 <- gopher(genes=T89.CTW4.CTW8$all, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW4.CTW8)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW4.CTW8$go,file = here("data/analysis/GO/T89.CTW4.CTW8-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW4.CTW8$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW4.CTW8-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.T89.CTW4.CTW8.up <- gopher(genes=T89.CTW4.CTW8$up, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW4.CTW8)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW4.CTW8.up$go,file = here("data/analysis/GO/T89.CTW4.CTW8.up-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW4.CTW8.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW4.CTW8.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.T89.CTW4.CTW8.dn <- gopher(genes=T89.CTW4.CTW8$dn, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW4.CTW8)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW4.CTW8.dn$go,file = here("data/analysis/GO/T89.CTW4.CTW8.dn-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW4.CTW8.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW4.CTW8.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment T89 CTW8 _vs._ LDD7
#' #' all DE genes
#' enr.T89.CTW8.LDD7 <- gopher(genes=T89.CTW8.LDD7$all, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW8.LDD7)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW8.LDD7$go,file = here("data/analysis/GO/T89.CTW8.LDD7-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW8.LDD7$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW8.LDD7-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.T89.CTW8.LDD7.up <- gopher(genes=T89.CTW8.LDD7$up, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW8.LDD7)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW8.LDD7.up$go,file = here("data/analysis/GO/T89.CTW8.LDD7.up-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW8.LDD7.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW8.LDD7.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.T89.CTW8.LDD7.dn <- gopher(genes=T89.CTW8.LDD7$dn, background=rownames(vst)[rowSums(counts(ddsEx.T89.CTW8.LDD7)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.T89.CTW8.LDD7.dn$go,file = here("data/analysis/GO/T89.CTW8.LDD7.dn-GO-enrichment.csv"))
#' write_delim(enr.T89.CTW8.LDD7.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/T89.CTW8.LDD7.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' 
#' 
#' #' ### GO enrichment FT1 LD _vs._ SDW15
#' #' all DE genes
#' enr.FT1.LD.SDW15 <- gopher(genes=FT1.LD.SDW15$all, background=rownames(vst)[rowSums(counts(ddsEx.FT1.LD.SDW15)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.LD.SDW15$go,file = here("data/analysis/GO/FT1.LD.SDW15-GO-enrichment.csv"))
#' write_delim(enr.FT1.LD.SDW15$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.LD.SDW15-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.FT1.LD.SDW15.up <- gopher(genes=FT1.LD.SDW15$up, background=rownames(vst)[rowSums(counts(ddsEx.FT1.LD.SDW15)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.LD.SDW15.up$go,file = here("data/analysis/GO/FT1.LD.SDW15.up-GO-enrichment.csv"))
#' write_delim(enr.FT1.LD.SDW15.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.LD.SDW15.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.FT1.LD.SDW15.dn <- gopher(genes=FT1.LD.SDW15$dn, background=rownames(vst)[rowSums(counts(ddsEx.FT1.LD.SDW15)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.LD.SDW15.dn$go,file = here("data/analysis/GO/FT1.LD.SDW15-dn-GO-enrichment.csv"))
#' write_delim(enr.FT1.LD.SDW15.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.LD.SDW15.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment FT1 SDW15 _vs._ CTW2
#' #' all DE genes
#' enr.FT1.SDW15.CTW2 <- gopher(genes=FT1.SDW15.CTW2$all, background=rownames(vst)[rowSums(counts(ddsEx.FT1.SDW15.CTW2)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.SDW15.CTW2$go,file = here("data/analysis/GO/FT1.SDW15.CTW2-GO-enrichment.csv"))
#' write_delim(enr.FT1.SDW15.CTW2$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.SDW15.CTW2-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.FT1.SDW15.CTW2.up <- gopher(genes=FT1.SDW15.CTW2$up, background=rownames(vst)[rowSums(counts(ddsEx.FT1.SDW15.CTW2)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.SDW15.CTW2.up$go,file = here("data/analysis/GO/FT1.SDW15.CTW2.up-GO-enrichment.csv"))
#' write_delim(enr.FT1.SDW15.CTW2.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.SDW15.CTW2.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.FT1.SDW15.CTW2.dn <- gopher(genes=FT1.SDW15.CTW2$dn, background=rownames(vst)[rowSums(counts(ddsEx.FT1.SDW15.CTW2)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.SDW15.CTW2.dn$go,file = here("data/analysis/GO/FT1.SDW15.CTW2.dn-GO-enrichment.csv"))
#' write_delim(enr.FT1.SDW15.CTW2.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.SDW15.CTW2.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment FT1 CTW2 _vs._ CTW4
#' #' all DE genes
#' enr.FT1.CTW2.CTW4 <- gopher(genes=FT1.CTW2.CTW4$all, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW2.CTW4)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW2.CTW4$go,file = here("data/analysis/GO/FT1.CTW2.CTW4-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW2.CTW4$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW2.CTW4-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.FT1.CTW2.CTW4.up <- gopher(genes=FT1.CTW2.CTW4$up, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW2.CTW4)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW2.CTW4.up$go,file = here("data/analysis/GO/FT1.CTW2.CTW4.up-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW2.CTW4.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW2.CTW4.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.FT1.CTW2.CTW4.dn <- gopher(genes=FT1.CTW2.CTW4$dn, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW2.CTW4)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW2.CTW4.dn$go,file = here("data/analysis/GO/FT1.CTW2.CTW4.dn-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW2.CTW4.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW2.CTW4.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment FT1 CTW4 _vs._ CTW8
#' #' all DE genes
#' enr.FT1.CTW4.CTW8 <- gopher(genes=FT1.CTW4.CTW8$all, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW4.CTW8)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW4.CTW8$go,file = here("data/analysis/GO/FT1.CTW4.CTW8-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW4.CTW8$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW4.CTW8-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.FT1.CTW4.CTW8.up <- gopher(genes=FT1.CTW4.CTW8$up, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW4.CTW8)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW4.CTW8.up$go,file = here("data/analysis/GO/FT1.CTW4.CTW8.up-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW4.CTW8.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW4.CTW8.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.FT1.CTW4.CTW8.dn <- gopher(genes=FT1.CTW4.CTW8$dn, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW4.CTW8)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW4.CTW8.dn$go,file = here("data/analysis/GO/FT1.CTW4.CTW8.dn-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW4.CTW8.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW4.CTW8.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment FT1 CTW8 _vs._ LDD7
#' #' all DE genes
#' enr.FT1.CTW8.LDD7 <- gopher(genes=FT1.CTW8.LDD7$all, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW8.LDD7)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW8.LDD7$go,file = here("data/analysis/GO/FT1.CTW8.LDD7-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW8.LDD7$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW8.LDD7-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.FT1.CTW8.LDD7.up <- gopher(genes=FT1.CTW8.LDD7$up, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW8.LDD7)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW8.LDD7.up$go,file = here("data/analysis/GO/FT1.CTW8.LDD7.up-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW8.LDD7.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW8.LDD7.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.FT1.CTW8.LDD7.dn <- gopher(genes=FT1.CTW8.LDD7$dn, background=rownames(vst)[rowSums(counts(ddsEx.FT1.CTW8.LDD7)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.FT1.CTW8.LDD7.dn$go,file = here("data/analysis/GO/FT1.CTW8.LDD7.dn-GO-enrichment.csv"))
#' write_delim(enr.FT1.CTW8.LDD7.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/FT1.CTW8.LDD7.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #'
#' #' ## GO enrichment for T89 vs FT1 at different time points
#' #' ### GO enrichment LD T89 _vs._ FT1
#' #' all DE genes
#' enr.LD.T89.FT1 <- gopher(genes=LD.T89.FT1$all, background=rownames(vst)[rowSums(counts(ddsEx.LD.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.LD.T89.FT1$go,file = here("data/analysis/GO/LD.T89.FT1-GO-enrichment.csv"))
#' write_delim(enr.LD.T89.FT1$go[,c("id","padj")],
#'             file = here("data/analysis/GO/LD.T89.FT1-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.LD.T89.FT1.up <- gopher(genes=LD.T89.FT1$up, background=rownames(vst)[rowSums(counts(ddsEx.LD.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.LD.T89.FT1.up$go,file = here("data/analysis/GO/LD.T89.FT1.up-GO-enrichment.csv"))
#' write_delim(enr.LD.T89.FT1.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/LD.T89.FT1.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.LD.T89.FT1.dn <- gopher(genes=LD.T89.FT1$dn, background=rownames(vst)[rowSums(counts(ddsEx.LD.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.LD.T89.FT1.dn$go,file = here("data/analysis/GO/LD.T89.FT1.dn-GO-enrichment.csv"))
#' write_delim(enr.LD.T89.FT1.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/LD.T89.FT1.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment SDW15 T89 _vs._ FT1
#' #' all DE genes
#' enr.SDW15.T89.FT1 <- gopher(genes=SDW15.T89.FT1$all, background=rownames(vst)[rowSums(counts(ddsEx.SDW15.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.SDW15.T89.FT1$go,file = here("data/analysis/GO/SDW15.T89.FT1-GO-enrichment.csv"))
#' write_delim(enr.SDW15.T89.FT1$go[,c("id","padj")],
#'             file = here("data/analysis/GO/SDW15.T89.FT1-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.SDW15.T89.FT1.up <- gopher(genes=SDW15.T89.FT1$up, background=rownames(vst)[rowSums(counts(ddsEx.SDW15.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.SDW15.T89.FT1.up$go,file = here("data/analysis/GO/SDW15.T89.FT1.up-GO-enrichment.csv"))
#' write_delim(enr.SDW15.T89.FT1.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/SDW15.T89.FT1.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.SDW15.T89.FT1.dn <- gopher(genes=SDW15.T89.FT1$dn, background=rownames(vst)[rowSums(counts(ddsEx.SDW15.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.SDW15.T89.FT1.dn$go,file = here("data/analysis/GO/SDW15.T89.FT1.dn-GO-enrichment.csv"))
#' write_delim(enr.SDW15.T89.FT1.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/SDW15.T89.FT1.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment CTW2 T89 _vs._ FT1
#' #' all DE genes
#' enr.CTW2.T89.FT1 <- gopher(genes=CTW2.T89.FT1$all, background=rownames(vst)[rowSums(counts(ddsEx.CTW2.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW2.T89.FT1$go,file = here("data/analysis/GO/CTW2.T89.FT1-GO-enrichment.csv"))
#' write_delim(enr.CTW2.T89.FT1$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW2.T89.FT1-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.CTW2.T89.FT1.up <- gopher(genes=CTW2.T89.FT1$up, background=rownames(vst)[rowSums(counts(ddsEx.CTW2.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW2.T89.FT1.up$go,file = here("data/analysis/GO/CTW2.T89.FT1.up-GO-enrichment.csv"))
#' write_delim(enr.CTW2.T89.FT1.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW2.T89.FT1.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.CTW2.T89.FT1.dn <- gopher(genes=CTW2.T89.FT1$dn, background=rownames(vst)[rowSums(counts(ddsEx.CTW2.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW2.T89.FT1.dn$go,file = here("data/analysis/GO/CTW2.T89.FT1.dn-GO-enrichment.csv"))
#' write_delim(enr.CTW2.T89.FT1.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW2.T89.FT1.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' #' ### GO enrichment CTW4 T89 _vs._ FT1
#' #' all DE genes
#' enr.CTW4.T89.FT1 <- gopher(genes=CTW4.T89.FT1$all, background=rownames(vst)[rowSums(counts(ddsEx.CTW4.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW4.T89.FT1$go,file = here("data/analysis/GO/CTW4.T89.FT1-GO-enrichment.csv"))
#' write_delim(enr.CTW4.T89.FT1$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW4.T89.FT1-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.CTW4.T89.FT1.up <- gopher(genes=CTW4.T89.FT1$up, background=rownames(vst)[rowSums(counts(ddsEx.CTW4.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW4.T89.FT1.up$go,file = here("data/analysis/GO/CTW4.T89.FT1.up-GO-enrichment.csv"))
#' write_delim(enr.CTW4.T89.FT1.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW4.T89.FT1.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.CTW4.T89.FT1.dn <- gopher(genes=CTW4.T89.FT1$dn, background=rownames(vst)[rowSums(counts(ddsEx.CTW4.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW4.T89.FT1.dn$go,file = here("data/analysis/GO/CTW4.T89.FT1.dn-GO-enrichment.csv"))
#' write_delim(enr.CTW4.T89.FT1.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW4.T89.FT1.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment CTW8 T89 _vs._ FT1
#' #' all DE genes
#' enr.CTW8.T89.FT1 <- gopher(genes=CTW8.T89.FT1$all, background=rownames(vst)[rowSums(counts(ddsEx.CTW8.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW8.T89.FT1$go,file = here("data/analysis/GO/CTW8.T89.FT1-GO-enrichment.csv"))
#' write_delim(enr.CTW8.T89.FT1$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW8.T89.FT1-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.CTW8.T89.FT1.up <- gopher(genes=CTW8.T89.FT1$up, background=rownames(vst)[rowSums(counts(ddsEx.CTW8.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW8.T89.FT1.up$go,file = here("data/analysis/GO/CTW8.T89.FT1.up-GO-enrichment.csv"))
#' write_delim(enr.CTW8.T89.FT1.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW8.T89.FT1.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.CTW8.T89.FT1.dn <- gopher(genes=CTW8.T89.FT1$dn, background=rownames(vst)[rowSums(counts(ddsEx.CTW8.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.CTW8.T89.FT1.dn$go,file = here("data/analysis/GO/CTW8.T89.FT1.dn-GO-enrichment.csv"))
#' write_delim(enr.CTW8.T89.FT1.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/CTW8.T89.FT1.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #'
#' #' ### GO enrichment LDD7 T89 _vs._ FT1
#' #' all DE genes
#' enr.LDD7.T89.FT1 <- gopher(genes=LDD7.T89.FT1$all, background=rownames(vst)[rowSums(counts(ddsEx.LDD7.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.LDD7.T89.FT1$go,file = here("data/analysis/GO/LDD7.T89.FT1-GO-enrichment.csv"))
#' write_delim(enr.LDD7.T89.FT1$go[,c("id","padj")],
#'             file = here("data/analysis/GO/LDD7.T89.FT1-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' up-regulated genes
#' enr.LDD7.T89.FT1.up <- gopher(genes=LDD7.T89.FT1$up, background=rownames(vst)[rowSums(counts(ddsEx.LDD7.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.LDD7.T89.FT1.up$go,file = here("data/analysis/GO/LDD7.T89.FT1.up-GO-enrichment.csv"))
#' write_delim(enr.LDD7.T89.FT1.up$go[,c("id","padj")],
#'             file = here("data/analysis/GO/LDD7.T89.FT1.up-for-REVIGO.txt"),
#'             col_names = FALSE)
#' #' down-regulated genes
#' enr.LDD7.T89.FT1.dn <- gopher(genes=LDD7.T89.FT1$dn, background=rownames(vst)[rowSums(counts(ddsEx.LDD7.T89.FT1)) != 0], task=c("go"), url="potra2")
#' write_csv(enr.LDD7.T89.FT1.dn$go,file = here("data/analysis/GO/LDD7.T89.FT1.dn-GO-enrichment.csv"))
#' write_delim(enr.LDD7.T89.FT1.dn$go[,c("id","padj")],
#'             file = here("data/analysis/GO/LDD7.T89.FT1.dn-for-REVIGO.txt"),
#'             col_names = FALSE)
#' 


#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```

