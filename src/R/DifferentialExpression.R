#' ---
#' title: "Differential Expression"
#' author: "Nicolas Delhomme & Alice Marcon"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup

#' * Libraries
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
    library(here)
    library(dplyr)
    library(magrittr)
    library(plotly)
    library(readr)
})

#' * Helper files
suppressMessages({
    source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    source(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
    message(paste("Plotting",gene_id))
    sel <- grepl(gene_id,rownames(vst))
    stopifnot(sum(sel)==1)

    p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                          data.frame(value=vst[sel,])),
                aes(x=Treatment,y=value,col=Genotype,
                    group=Genotype,
                    text=NGI_ID)) +
        geom_point() + geom_smooth() +
        scale_y_continuous(name="VST expression") + 
        ggtitle(label=paste("Expression for: ",gene_id))
    
    suppressMessages(suppressWarnings(plot(p)))
    return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds,vst,contrast,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),...){
    
    # get the filter
    if(!is.null(match.arg(filter))){
        filter <- rowMedians(counts(dds,normalized=TRUE))
        message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    # validation
    if(length(contrast)==1){
        res <- results(dds,name=contrast,filter = filter,lfcThreshold=lfc,alpha=padj)
    } else {
        res <- results(dds,contrast=contrast,filter = filter,lfcThreshold=lfc,alpha=padj)
    }
    
    stopifnot(length(sample_sel)==ncol(vst))
    
    if(plot){
        par(mar=c(5,5,5,5))
        volcanoPlot(res)
        par(mar=mar)
    }
    
    # a look at independent filtering
    if(plot){
        plot(metadata(res)$filterNumRej,
             type="b", ylab="number of rejections",
             xlab="quantiles of filter")
        lines(metadata(res)$lo.fit, col="red")
        abline(v=metadata(res)$filterTheta)
    }
    
    if(verbose){
        message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                        round(metadata(res)$filterThreshold,digits=5),
                        names(metadata(res)$filterThreshold)))
        
        max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
        message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                        round(max.theta*100,digits=5),
                        round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
    }
    
    if(plot){
        qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
        dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                          qtl.exp=qtl.exp,
                          number.degs=sapply(lapply(qtl.exp,function(qe){
                              res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                                  ! is.na(res$padj) & res$baseMean >= qe
                          }),sum))
        if(debug){
            plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("base mean expression") +
                     geom_hline(yintercept=expression_cutoff,
                                linetype="dotted",col="red"))
        
            p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
                geom_line() + geom_point() +
                scale_x_continuous("quantiles of expression") + 
                scale_y_log10("base mean expression") + 
                geom_hline(yintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
                     geom_line() + geom_point() +
                     geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes"))
            
            plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Cumulative number of DE genes"))
            
            plot(ggplot(data.frame(x=dat$thetas[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("quantiles of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                     geom_line() + geom_point() +
                     scale_x_continuous("base mean of expression") + 
                     scale_y_continuous("Number of DE genes per interval"))
            
            p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                                   y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
                geom_line() + geom_point() +
                scale_x_log10("base mean of expression") + 
                scale_y_continuous("Number of DE genes per interval") + 
                geom_vline(xintercept=expression_cutoff,
                           linetype="dotted",col="red")
            suppressMessages(suppressWarnings(plot(p)))
        }
    }
    
    sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
        res$baseMean >= expression_cutoff
    
    if(verbose){
      message(sprintf(paste(
        ifelse(sum(sel)==1,
               "There is %s gene that is DE",
               "There are %s genes that are DE"),
        "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
        sum(sel),padj,
        lfc,expression_cutoff))
    }
    
    # proceed only if there are DE genes
    if(sum(sel) > 0){
        val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
        if (sum(val) >0){
          warning(sprintf(paste(
            ifelse(sum(val)==1,
                   "There is %s DE gene that has",
                   "There are %s DE genes that have"),
            "no vst expression in the selected samples"),sum(val)))
          sel[sel][val] <- FALSE
        } 

        if(export){
            if(!dir.exists(default_dir)){
                dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
            }
            write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
            write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
        }
        if(plot & sum(sel)>1){
            heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                      distfun = pearson.dist,
                      hclustfun = function(X){hclust(X,method="ward.D2")},
                      trace="none",col=hpal,labRow = FALSE,
                      labCol=labels[sample_sel],...
            )
        }
    }
    return(list(all=rownames(res[sel,]),
                up=rownames(res[sel & res$log2FoldChange > 0,]),
                dn=rownames(res[sel & res$log2FoldChange < 0,])))
}

#' 3. extract and plot the enrichment results
extractEnrichmentResults <- function(enrichment,task="go",
                                     diff.exp=c("all","up","dn"),
                                     go.namespace=c("BP","CC","MF"),
                                     genes=NULL,export=TRUE,plot=TRUE,
                                     default_dir=here("data/analysis/DE"),
                                     default_prefix="DE",
                                     url="athaliana"){
    # process args
    diff.exp <- match.arg(diff.exp)
    de <- ifelse(diff.exp=="all","none",
                 ifelse(diff.exp=="dn","down",diff.exp))

    # sanity
    if( is.null(enrichment[[task]]) | length(enrichment[[task]]) == 0){
        message(paste("No enrichment for",task))
    } else {

        # write out
        if(export){
            write_tsv(enrichment[[task]],
                      file=here(default_dir,
                                paste0(default_prefix,"-genes_GO-enrichment.tsv")))
            if(!is.null(genes)){
                write_tsv(
                    enrichedTermToGenes(genes=genes,terms=enrichment[[task]]$id,url=url,mc.cores=16L),
                    file=here(default_dir,
                              paste0(default_prefix,"-enriched-term-to-genes.tsv"))
                )
            }
        }
        
        if(plot){
            sapply(go.namespace,function(ns){
                titles <- c(BP="Biological Process",
                            CC="Cellular Component",
                            MF="Molecular Function")
                suppressWarnings(tryCatch({plotEnrichedTreemap(enrichment,enrichment=task,
                                                               namespace=ns,
                                                               de=de,title=paste(default_prefix,titles[ns]))},
                                          error = function(e) {
                                              message(paste("Treemap plot failed for",ns, 
                                                            "because of:",e))
                                              return(NULL)
                                          }))
            })
        }
    }
}

#' * Data
load("/mnt/picea/home/dandre/Git/aspen-FTL1-growth-cessation/analysis/FT1_CRISPR_V2_app/dds_corrected.rda")

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
                                         "FT1_LDD7_1", "FT1_LDD7_2", "FT1_LDD7_3", "FT1_LDD7_4", "FT1_LDD7_5", "FT1_LDD7_6"))

dds$Genotype <- droplevels(dds$Genotype)
dds$Batch <- factor(substr(dds$NGI_ID,8,8))

#' Finally update the design
design(dds) <- ~Genotype * Treatment

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Save dds and vst
save(vst,file=here("../FT1 CRISPR/vst_corrected.rda"))
save(dds,file=here("../FT1 CRISPR/dds_corrected.rda"))


#' ## Gene of interest
#' remove line 736 and outliers
sel.T89.FT1 <- dds$Genotype != "736" & ! dds$NGI_ID %in% c("P12108_101", "P12108_221", "P12108_231")

#' Heatmap of "all" genes
#' Create new order (by time point)
sel2 <- c()
sel2 <- order(dds$UserID)

#re-order the columns of the dds
dds <- dds[,sel2]
# and the vst
vst <- vst[,sel2]
# double checking
stopifnot(colnames(vst)==colnames(dds))

save(dds,vst,file=here("../FT1 CRISPR/dds_vst_in_order_for_heatmap.rda"))

## Heatmap of all genes 
#' Same as in biological QA, but only with T89 and FT1 and without outliers
conds <- factor(paste(dds$Genotype,dds$Treatment)[sel.T89.FT1],
                levels=c("T89 LD", "FT1 LD" ,
                         "T89 SDW15","FT1 SDW15",
                         "T89 CTW2", "FT1 CTW2",
                         "T89 CTW4","FT1 CTW4",   
                         "T89 CTW8","FT1 CTW8",
                         "T89 LDD7",    "FT1 LDD7"))

sels <- rangeFeatureSelect(counts=vst[,sel.T89.FT1],
                           conditions=conds,
                           nrep=3)

hm <- heatmap.2(t(scale(t(vst[sels[[6]],sel.T89.FT1]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,cexCol=.6,
                col=hpal)

#select in the vst only the rows for few genes 
genes_for_heatmap <- read_table2("~/FT1 CRISPR/genes_for_heatmap.txt")
genes <- genes_for_heatmap$ID

# dataframe with the name of the labels for the heatmap
genes_for_heatmap_p <- read_csv("~/FT1 CRISPR/genes_for_heatmap_p.txt")

# create the vst for the genes that have been selected
vst_small <- vst[genes,]
vst_small <- vst_small[rowSums(vst_small) > 0,]

separatorrow<- c(5,8,25,28,31,39)
separatorcol <- c(11,23,34,45,57)

hm <- heatmap.2(t(scale(t(vst_small[,sel.T89.FT1]))),
                Colv=FALSE, dendrogram="none",
                Rowv = FALSE,
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow =genes_for_heatmap_p$`Gene ID`, 
                trace = "none", cexRow = .6,
                labCol = conds,cexCol=.6,
                margins = c(5,9),
                colsep = separatorcol,
                rowsep = separatorrow,
                sepcolor = "black",
                sepwidth = c(0.02,0.02),
                col=hpal,key = FALSE,keysize = .2)


# expression is low overall
# hm <- heatmap.2(vst_small[,sel.T89.FT1],
#                 Colv=FALSE, dendrogram="row",
#                 labRow = genes_for_heatmap$Gene,trace = "none", cexRow = .4,
#                 labCol = conds,cexCol=.6,
#                 colsep = separatorcol,
#                 rowsep = separatorrow,
#                 sepcolor = "white",
#                 sepwidth = c(0.05,0.05),
#                 col=hpal)
# 
# cutree(hm$rowDendrogram,h=3)


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

#' # Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)

#' How many genes expressed at different timepoint
res <- results(dds)
#' Percent of genes non informative (noise) and corresponding BaseMean cutoff
#' (for all samples)
metadata(res)$filterThreshold
sel <- !is.na(res$padj)
#' For T89 and timepoint LD
mean(colSums(vst[sel,dds$Treatment=="LD" & dds$Genotype== "T89"] > 0))


#'DE analysis between ft1 and T89 at each timepoint

#'ft1 vs T89 at LD
FT1_T89_LD <- extract_results(dds=dds,
                                vst=vst,
                                contrast="Genotype_FT1_vs_T89",
                                export = FALSE,
                                sample_sel = dds$Treatment == "LD",
                                labels = dds$UserID,
                                cexCol=0.7)
#'ft1 vs T89 at SDW15
contrast.sel <- rep(0,length(resultsNames(dds)))
names(contrast.sel) <- resultsNames(dds)
contrast.sel["Genotype_FT1_vs_T89"] <- 1
contrast.sel["GenotypeFT1.TreatmentSDW15"] <- 1
FT1_T89_SDW15 <- extract_results(dds=dds,
                                vst=vst,
                                contrast=contrast.sel,
                                export = FALSE,
                                sample_sel = dds$Treatment == "SDW15",
                                labels = dds$UserID,
                                cexCol=0.7)

#'ft1 vs T89 at CTW2
contrast.sel["GenotypeFT1.TreatmentSDW15"] <- 0
contrast.sel["GenotypeFT1.TreatmentCTW2"] <- 1
FT1_T89_CTW2 <- extract_results(dds=dds,
                                 vst=vst,
                                 contrast=contrast.sel,
                                 export = FALSE,
                                 sample_sel = dds$Treatment == "CTW2",
                                 labels = dds$UserID,
                                 cexCol=0.7)

#'ft1 vs T89 at CTW4
contrast.sel["GenotypeFT1.TreatmentCTW2"] <- 0
contrast.sel["GenotypeFT1.TreatmentCTW4"] <- 1
FT1_T89_CTW4 <- extract_results(dds=dds,
                                vst=vst,
                                contrast=contrast.sel,
                                export = FALSE,
                                sample_sel = dds$Treatment == "CTW4",
                                labels = dds$UserID,
                                cexCol=0.7)

#'ft1 vs T89 at CTW8
contrast.sel["GenotypeFT1.TreatmentCTW4"] <- 0
contrast.sel["GenotypeFT1.TreatmentCTW8"] <- 1
FT1_T89_CTW8 <- extract_results(dds=dds,
                                vst=vst,
                                contrast=contrast.sel,
                                export = FALSE,
                                sample_sel = dds$Treatment == "CTW8",
                                labels = dds$UserID,
                                cexCol=0.7)
#'ft1 vs T89 at LDD7
contrast.sel["GenotypeFT1.TreatmentCTW8"] <- 0
contrast.sel["GenotypeFT1.TreatmentLDD7"] <- 1
FT1_T89_LDD7 <- extract_results(dds=dds,
                                vst=vst,
                                contrast=contrast.sel,
                                export = FALSE,
                                sample_sel = dds$Treatment == "LDD7",
                                labels = dds$UserID,
                                cexCol=0.7)


#' Change to have CTW4 as the baseline 
dds$Treatment <- relevel(dds$Treatment,"CTW4")
#' # Re-run Differential Expression
dds <- DESeq(dds)
#' Check the different contrasts
resultsNames(dds)

#' Time effect T89 CTW8 _vs._ CTW4
T89_CTW8_CTW4 <- extract_results(dds=dds,
                                vst=vst,
                                contrast="Treatment_CTW8_vs_CTW4",
                                export = FALSE,
                                sample_sel = dds$Treatment %in% c("CTW4","CTW8") & 
                                  dds$Genotype == "T89",
                                labels = dds$UserID,
                                cexCol=0.7)

#' Time effect FT1 CTW8 _vs._ CTW4
contrast.sel <- rep(0,length(resultsNames(dds)))
names(contrast.sel) <- resultsNames(dds)
contrast.sel["Treatment_CTW8_vs_CTW4"] <- 1
contrast.sel["GenotypeFT1.TreatmentCTW8"] <- 1

FT1_CTW8_CTW4 <- extract_results(dds=dds,
                                 vst=vst,
                                 contrast=contrast.sel,
                                 export = FALSE,
                                 sample_sel = dds$Treatment %in% c("CTW4","CTW8") & 
                                   dds$Genotype == "FT1",
                                 labels = dds$UserID,
                                 cexCol=0.7)

#' ### DE contrast for LDD7
#' Change to have CTW8 as the baseline 
#' dds$Treatment <- relevel(dds$Treatment,"CTW8")
#' #' # Re-run Differential Expression
#' dds <- DESeq(dds)
#' #' Check the different contrasts
#' resultsNames(dds)

#' #' Time effect T89 LDD7 _vs._ CTW8
#' T89_LDD7_CTW8 <- extract_results(dds=dds,
#'                                  vst=vst,
#'                                  contrast="Treatment_LDD7_vs_CTW8",
#'                                  export = FALSE,
#'                                  sample_sel = dds$Treatment %in% c("CTW8","LDD7") & 
#'                                    dds$Genotype == "T89",
#'                                  labels = dds$UserID,
#'                                  cexCol=0.7)
#' 
#' #' Time effect FT1 LDD7 _vs._ CTW8
#' contrast.sel <- rep(0,length(resultsNames(dds)))
#' names(contrast.sel) <- resultsNames(dds)
#' contrast.sel["Treatment_LDD7_vs_CTW8"] <- 1
#' contrast.sel["GenotypeFT1.TreatmentLDD7"] <- 1
#' FT1_LDD7_CTW8 <- extract_results(dds=dds,
#'                                  vst=vst,
#'                                  contrast=contrast.sel,
#'                                  export = FALSE,
#'                                  sample_sel = dds$Treatment %in% c("CTW8","LDD7") & 
#'                                    dds$Genotype == "FT1",
#'                                  labels = dds$UserID,
#'                                  cexCol=0.7)


#' ### Venn Diagram choose the contrast 
res.list <- list(FT1=FT1_CTW8_CTW4,T89=T89_CTW8_CTW4)

#' #### All DE genes for ft1 and T89 between CTW8 and CTW4
grid.newpage()
grid.draw(venn.diagram(lapply(res.list,"[[","all"),
                       NULL,
                       fill=pal[1:2]))

## GO enrichment for DE genes at CTW8 between FT1 and T89
#' #' all DE genes
enr.CTW8.T89.FT1 <- gopher(genes=FT1_T89_CTW8$all, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.CTW8.T89.FT1$go,file = here("/mnt/picea/home/amarcon/FT1 CRISPR/GO/CTW8.T89.FT1-GO-enrichment.csv"))
write_delim(enr.CTW8.T89.FT1$go[,c("id","padj")],
            file = here("/mnt/picea/home/amarcon/FT1 CRISPR/GO/CTW8.T89.FT1-for-REVIGO.txt"),
            col_names = FALSE)
#' #' up-regulated genes
enr.CTW8.T89.FT1.up <- gopher(genes=FT1_T89_CTW8$up, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.CTW8.T89.FT1.up$go,file = here("/mnt/picea/home/amarcon/FT1 CRISPR/GO/CTW8.T89.FT1.up-GO-enrichment.csv"))
write_delim(enr.CTW8.T89.FT1.up$go[,c("id","padj")],
           file = here("/mnt/picea/home/amarcon/FT1 CRISPR/GO/CTW8.T89.FT1.up-for-REVIGO.txt"),
           col_names = FALSE)
#' #' down-regulated genes
enr.CTW8.T89.FT1.dn <- gopher(genes=FT1_T89_CTW8$dn, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.CTW8.T89.FT1.dn$go,file = here("/mnt/picea/home/amarcon/FT1 CRISPR/GO/CTW8.T89.FT1.dn-GO-enrichment.csv"))
write_delim(enr.CTW8.T89.FT1.dn$go[,c("id","padj")],
           file = here("/mnt/picea/home/amarcon/FT1 CRISPR/GO/CTW8.T89.FT1.dn-for-REVIGO.txt"),
          col_names = FALSE)

## GO enrichment for DE genes in T89 between CTW8 and CTW4

enr.T89.CTW8.CTW4 <- gopher(genes=T89_CTW8_CTW4$all, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.T89.CTW8.CTW4$go,file = here("../FT1 CRISPR/GO/T89.CTW8.CTW4-GO-enrichment.csv"))
write_delim(enr.T89.CTW8.CTW4$go[,c("id","padj")],
            file = here("../FT1 CRISPR/GO/T89.CTW8.CTW4-for-REVIGO.txt"),
            col_names = FALSE)

enr.T89.CTW8.CTW4.up <- gopher(genes=T89_CTW8_CTW4$up, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.T89.CTW8.CTW4.up$go,file = here("../FT1 CRISPR/GO/T89.CTW8.CTW4.up-GO-enrichment.csv"))
write_delim(enr.T89.CTW8.CTW4.up$go[,c("id","padj")],
            file = here("../FT1 CRISPR/GO/T89.CTW8.CTW4.up-for-REVIGO.txt"),
            col_names = FALSE)

enr.T89.CTW8.CTW4.dn <- gopher(genes=T89_CTW8_CTW4$dn, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.T89.CTW8.CTW4.dn$go,file = here("../FT1 CRISPR/GO/T89.CTW8.CTW4.dn-GO-enrichment.csv"))
write_delim(enr.T89.CTW8.CTW4.dn$go[,c("id","padj")],
            file = here("../FT1 CRISPR/GO/T89.CTW8.CTW4.dn-for-REVIGO.txt"),
            col_names = FALSE)

## GO enrichment for DE genes in FT1 between CTW8 and CTW4

enr.FT1.CTW8.CTW4 <- gopher(genes=FT1_CTW8_CTW4$all, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.FT1.CTW8.CTW4$go,file = here("../FT1 CRISPR/GO/FT1.CTW8.CTW4-GO-enrichment.csv"))
write_delim(enr.FT1.CTW8.CTW4$go[,c("id","padj")],
            file = here("../FT1 CRISPR/GO/FT1.CTW8.CTW4-for-REVIGO.txt"),
            col_names = FALSE)

enr.FT1.CTW8.CTW4.up <- gopher(genes=FT1_CTW8_CTW4$up, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.FT1.CTW8.CTW4.up$go,file = here("../FT1 CRISPR/GO/FT1.CTW8.CTW4.up-GO-enrichment.csv"))
write_delim(enr.FT1.CTW8.CTW4.up$go[,c("id","padj")],
            file = here("../FT1 CRISPR/GO/FT1.CTW8.CTW4.up-for-REVIGO.txt"),
            col_names = FALSE)

enr.FT1.CTW8.CTW4.dn <- gopher(genes=FT1_CTW8_CTW4$dn, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
write_csv(enr.FT1.CTW8.CTW4.dn$go,file = here("../FT1 CRISPR/GO/FT1.CTW8.CTW4.dn-GO-enrichment.csv"))
write_delim(enr.FT1.CTW8.CTW4.dn$go[,c("id","padj")],
            file = here("../FT1 CRISPR/GO/FT1.CTW8.CTW4.dn-for-REVIGO.txt"),
            col_names = FALSE)


## GO enrichment for DE genes in FT1 between LDD7 and CTW8
# enr.FT1.LDD7.CTW8 <- gopher(genes=FT1_LDD7_CTW8$all, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.FT1.LDD7.CTW8$go,file = here("../FT1 CRISPR/GO/FT1.LDD7.CTW8-GO-enrichment.csv"))
# write_delim(enr.FT1.LDD7.CTW8$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/FT1.LDD7.CTW8-for-REVIGO.txt"),
#             col_names = FALSE)
# 
# enr.FT1.LDD7.CTW8up <- gopher(genes=FT1_LDD7_CTW8$up, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.FT1.LDD7.CTW8up$go,file = here("../FT1 CRISPR/GO/FT1.LDD7.CTW8.up-GO-enrichment.csv"))
# write_delim(enr.FT1.LDD7.CTW8up$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/FT1.LDD7.CTW8.up-for-REVIGO.txt"),
#             col_names = FALSE)
# 
# enr.FT1.LDD7.CTW8dn <- gopher(genes=FT1_LDD7_CTW8$dn, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.FT1.LDD7.CTW8dn$go,file = here("../FT1 CRISPR/GO/FT1.LDD7.CTW8.dn-GO-enrichment.csv"))
# write_delim(enr.FT1.LDD7.CTW8dn$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/FT1.LDD7.CTW8.dn-for-REVIGO.txt"),
#             col_names = FALSE)
## GO enrichment for DE genes in T89 between LDD7 and CTW8
# enr.T89.LDD7.CTW8 <- gopher(genes=T89_LDD7_CTW8$all, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.T89.LDD7.CTW8$go,file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8-GO-enrichment.csv"))
# write_delim(enr.T89.LDD7.CTW8$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8-for-REVIGO.txt"),
#             col_names = FALSE)
# 
# enr.T89.LDD7.CTW8up <- gopher(genes=T89_LDD7_CTW8$up, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# # write_csv(enr.T89.LDD7.CTW8up$go,file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.up-GO-enrichment.csv"))
# # write_delim(enr.T89.LDD7.CTW8up$go[,c("id","padj")],
# #             file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.up-for-REVIGO.txt"),
# #             col_names = FALSE)
# 
# enr.T89.LDD7.CTW8dn <- gopher(genes=T89_LDD7_CTW8$dn, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# # write_csv(enr.T89.LDD7.CTW8dn$go,file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.dn-GO-enrichment.csv"))
# # write_delim(enr.T89.LDD7.CTW8dn$go[,c("id","padj")],
# #             file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.dn-for-REVIGO.txt"),
# #             col_names = FALSE)

## GO enrichment for DE genes at LDD7 between FT1 and T89
# enr.FT1.T89.LDD7 <- gopher(genes=FT1_T89_LDD7$all, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.T89.LDD7.CTW8$go,file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8-GO-enrichment.csv"))
# write_delim(enr.T89.LDD7.CTW8$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8-for-REVIGO.txt"),
#             col_names = FALSE)
# 
# enr.FT1.T89.LDD7up <- gopher(genes=FT1_T89_LDD7$up, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.T89.LDD7.CTW8up$go,file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.up-GO-enrichment.csv"))
# write_delim(enr.T89.LDD7.CTW8up$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.up-for-REVIGO.txt"),
#             col_names = FALSE)
# 
# enr.FT1.T89.LDD7dn <- gopher(genes=FT1_T89_LDD7$dn, background=rownames(vst)[rowSums(counts(dds)) != 0], task=c("go"), url="potra2")
# write_csv(enr.T89.LDD7.CTW8dn$go,file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.dn-GO-enrichment.csv"))
# write_delim(enr.T89.LDD7.CTW8dn$go[,c("id","padj")],
#             file = here("../FT1 CRISPR/GO/T89.LDD7.CTW8.dn-for-REVIGO.txt"),
#             col_names = FALSE)
# extract go terms for LDD7 all  FT1 vs T89
# extractEnrichmentResults(enr.FT1.T89.LDD7,
#                          diff.exp="all",
#                          genes=FT1_T89_LDD7$all,
#                          default_dir=here("../FT1 CRISPR"),
#                          default_prefix="FT1.T89.LDD7.all",
#                          url="potra2")
#table(FT1_T89_LDD7_all_enriched_term_to_genes$GOID)
#setdiff(enr.FT1.LDD7.CTW8dn$go$name, enr.T89.LDD7.CTW8dn$go$name)
#intersect(enr.FT1.LDD7.CTW8dn$go$name, enr.T89.LDD7.CTW8dn$go$name)

## extract go terms for FT1 vs T89 at CTW8
# all
extractEnrichmentResults(enr.CTW8.T89.FT1,
                         diff.exp="all",
                         genes=FT1_T89_CTW8$all,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="CTW8.T89.FT1.all",
                         url="potra2")
#up
extractEnrichmentResults(enr.CTW8.T89.FT1.up,
                         diff.exp="up",
                         genes=FT1_T89_CTW8$up,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="CTW8.T89.FT1.up",
                         url="potra2")
#down
extractEnrichmentResults(enr.CTW8.T89.FT1.dn,
                         diff.exp="dn",
                         genes=FT1_T89_CTW8$dn,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="CTW8.T89.FT1.dn",
                         url="potra2")

## extract go terms for FT1 between CTW8 ad CTW4
# all
extractEnrichmentResults(enr.FT1.CTW8.CTW4,
                         diff.exp="all",
                         genes=FT1_CTW8_CTW4$all,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="FT1.CTW8.CTW4.all",
                         url="potra2")
#up
extractEnrichmentResults(enr.FT1.CTW8.CTW4.up,
                         diff.exp="up",
                         genes=FT1_CTW8_CTW4$up,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="FT1.CTW8.CTW4.up",
                         url="potra2")
#down
extractEnrichmentResults(enr.FT1.CTW8.CTW4.dn,
                         diff.exp="dn",
                         genes=FT1_CTW8_CTW4$dn,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="FT1.CTW8.CTW4.dn",
                         url="potra2")

## extract go terms for T89 between CTW8 ad CTW4
# all
extractEnrichmentResults(enr.T89.CTW8.CTW4,
                         diff.exp="all",
                         genes=T89_CTW8_CTW4$all,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="T89.CTW8.CTW4.all",
                         url="potra2")
#up
extractEnrichmentResults(enr.T89.CTW8.CTW4.up,
                         diff.exp="up",
                         genes=T89_CTW8_CTW4$up,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="T89.CTW8.CTW4.up",
                         url="potra2")
#down
extractEnrichmentResults(enr.T89.CTW8.CTW4.dn,
                         diff.exp="dn",
                         genes=T89_CTW8_CTW4$dn,
                         default_dir=here("../FT1 CRISPR"),
                         default_prefix="T89.CTW8.CTW4.dn",
                         url="potra2")


## annotation for all, up ,down regulated genes 

# annotation, keep only 6 columns and reduce to the gene level
annot <- read_tsv(here("reference/annotation/blast2go/Potra22_blast2go_GO_export.txt.gz"), 
                  show_col_types = FALSE  ) %>%  
  mutate(GeneID=sub("\\.[0-9]+$","",`Sequence Name`)) %>% 
  dplyr::select(GeneID,`Sequence Description`,`Enzyme Name`,`InterPro Name`,`Annotation GO ID`,`Annotation GO Term`) %>%
  filter(!duplicated(GeneID))


#' create a file txt with list of up and down regulated DE genes 
# up genes at CTW8 between FT1 vs T89
write.table(FT1_T89_CTW8$up, file = "CTW8.DEup.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#create dataframe with list of genes up annotated
up<- left_join(read_delim(here("../FT1 CRISPR/CTW8.DEup.txt"),
                            col_names = c("ID","GeneID"),skip = 1,delim=" ",
                            col_select = 2,show_col_types = FALSE),
                 annot,by="GeneID")
colnames(up) <- c("GENEID", "Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")

#add the expression for FT1 and T89 and the delta 
load(here("../FT1 CRISPR/vst_corrected.rda"))
load(here("../FT1 CRISPR/dds_corrected.rda"))

col.sel <- dds$Treatment=="CTW8"
vst <- tibble(as.data.frame(vst) %>%
                rownames_to_column("GENEID")) %>% 
  dplyr::select("GENEID",dds$NGI_ID[col.sel]) %>% 
  pivot_longer(cols = -GENEID,
               names_to = "Samples",
               values_to = "VST") %>% 
  mutate(Genotype=dds$Genotype[match(Samples,dds$NGI_ID)]) %>% 
  group_by(Genotype,GENEID) %>% summarise(VSTmedian = median(VST))

vst_FT1_CTW8 <- vst[vst$Genotype == "FT1",]
vst_T89_CTW8 <- vst[vst$Genotype == "T89",]

up<-  merge(vst_FT1_CTW8, up, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
up <- up[,c(1,3:8)]
colnames(up) <- c("GENEID","VSTmedian_FT1","Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
up<-  merge(vst_T89_CTW8, up, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
up <- up[,c(1,3:9)]
colnames(up) <- c("GENEID","VSTmedian_T89","VSTmedian_FT1","Sequence Description", "Enzyme Name", "InterPro Name","Annotation GO ID","Annotation GO Term")
up %<>% mutate(Delta=VSTmedian_FT1 - VSTmedian_T89)
up <- up[,c(1,4:8,2,3,9)]

write_tsv(up,file =  here("../FT1 CRISPR/CTW8.DEup.annotated.tsv"))

# down genes at CTW8 between FT1 vs T89
write.table(FT1_T89_CTW8$dn, file = "CTW8.DEdn.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#create dataframe with list of genes down annotated
down <- left_join(read_delim(here("../FT1 CRISPR/CTW8.DEdn.txt"),
                          col_names = c("ID","GeneID"),skip = 1,delim=" ",
                          col_select = 2,show_col_types = FALSE),
               annot,by="GeneID")
colnames(down) <- c("GENEID", "Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
down<-  merge(vst_FT1_CTW8, down, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
down <- down[,c(1,3:8)]
colnames(down) <- c("GENEID","VSTmedian_FT1","Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
down<-  merge(vst_T89_CTW8, down, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
down <- down[,c(1,3:9)]
colnames(down) <- c("GENEID","VSTmedian_T89","VSTmedian_FT1","Sequence Description", "Enzyme Name", "InterPro Name","Annotation GO ID","Annotation GO Term")
down %<>% mutate(Delta=VSTmedian_FT1 - VSTmedian_T89)
down <- down[,c(1,4:8,2,3,9)]

write_tsv(down,file =  here("../FT1 CRISPR/T89.DEdn.annotated.tsv"))


# up genes for FT1 between CTW8 and CTW4
write.table(FT1_CTW8_CTW4$up, file = "FT1.DEup.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
#create dataframe with list of genes up annotated
up<- left_join(read_delim(here("../FT1 CRISPR/FT1.DEup.txt"),
                          col_names = c("ID","GeneID"),skip = 1,delim=" ",
                          col_select = 2,show_col_types = FALSE),
               annot,by="GeneID")
colnames(up) <- c("GENEID", "Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")

#add the expression for FT1 CTW8 and CTW4
load(here("../FT1 CRISPR/vst_corrected.rda"))
load(here("../FT1 CRISPR/dds_corrected.rda"))

col.sel <- dds$Treatment=="CTW4"
vst <- tibble(as.data.frame(vst) %>%
                rownames_to_column("GENEID")) %>% 
  dplyr::select("GENEID",dds$NGI_ID[col.sel]) %>% 
  pivot_longer(cols = -GENEID,
               names_to = "Samples",
               values_to = "VST") %>% 
  mutate(Genotype=dds$Genotype[match(Samples,dds$NGI_ID)]) %>% 
  group_by(Genotype,GENEID) %>% summarise(VSTmedian = median(VST))

vst_FT1_CTW4 <- vst[vst$Genotype == "FT1",]

up<-  merge(vst_FT1_CTW8, up, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
up <- up[,c(1,3:8)]
colnames(up) <- c("GENEID","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
up<-  merge(vst_FT1_CTW4, up, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
up <- up[,c(1,3:9)]
colnames(up) <- c("GENEID","VSTmedian_CTW4","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name","Annotation GO ID","Annotation GO Term")
up %<>% mutate(Delta=VSTmedian_CTW8 - VSTmedian_CTW4)
up <- up[,c(1,4:8,2,3,9)]

write_tsv(up,file =  here("../FT1 CRISPR/FT1.DEup.annotated.tsv"))

# down for FT1 between CTW8 and CTW4
write.table(FT1_CTW8_CTW4$dn, file = "FT1.DEdn.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#create dataframe with list of genes down annotated
down <- left_join(read_delim(here("../FT1 CRISPR/FT1.DEdn.txt"),
                             col_names = c("ID","GeneID"),skip = 1,delim=" ",
                             col_select = 2,show_col_types = FALSE),
                  annot,by="GeneID")
colnames(down) <- c("GENEID", "Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
down<-  merge(vst_FT1_CTW8, down, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
down <- down[,c(1,3:8)]
colnames(down) <- c("GENEID","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
down<-  merge(vst_FT1_CTW4, down, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
down <- down[,c(1,3:9)]
colnames(down) <- c("GENEID","VSTmedian_CTW4","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name","Annotation GO ID","Annotation GO Term")
down %<>% mutate(Delta=VSTmedian_CTW8 - VSTmedian_CTW4)
down <- down[,c(1,4:8,2,3,9)]

write_tsv(down,file =  here("../FT1 CRISPR/FT1.DEdn.annotated.tsv"))



# up genes for T89 between CTW8 and CTW4
write.table(T89_CTW8_CTW4$up, file = "T89.DEup.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)
#create dataframe with list of genes up annotated
up<- left_join(read_delim(here("../FT1 CRISPR/T89.DEup.txt"),
                          col_names = c("ID","GeneID"),skip = 1,delim=" ",
                          col_select = 2,show_col_types = FALSE),
               annot,by="GeneID")
colnames(up) <- c("GENEID", "Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")

#add the expression for FT1 CTW8 and CTW4
vst_T89_CTW4 <- vst[vst$Genotype == "T89",]

up<-  merge(vst_T89_CTW8, up, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
up <- up[,c(1,3:8)]
colnames(up) <- c("GENEID","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
up<-  merge(vst_T89_CTW4, up, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
up <- up[,c(1,3:9)]
colnames(up) <- c("GENEID","VSTmedian_CTW4","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name","Annotation GO ID","Annotation GO Term")
up %<>% mutate(Delta=VSTmedian_CTW8 - VSTmedian_CTW4)
up <- up[,c(1,4:8,2,3,9)]

write_tsv(up,file =  here("../FT1 CRISPR/T89.DEup.annotated.tsv"))

# down for T89 between CTW8 and CTW4
write.table(T89_CTW8_CTW4$dn, file = "T89.DEdn.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

#create dataframe with list of genes down annotated
down <- left_join(read_delim(here("../FT1 CRISPR/T89.DEdn.txt"),
                             col_names = c("ID","GeneID"),skip = 1,delim=" ",
                             col_select = 2,show_col_types = FALSE),
                  annot,by="GeneID")
colnames(down) <- c("GENEID", "Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
down<-  merge(vst_T89_CTW8, down, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
down <- down[,c(1,3:8)]
colnames(down) <- c("GENEID","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name", "Annotation GO ID","Annotation GO Term")
down<-  merge(vst_T89_CTW4, down, by.y= "GENEID", all.x= FALSE, all.y= FALSE)
down <- down[,c(1,3:9)]
colnames(down) <- c("GENEID","VSTmedian_CTW4","VSTmedian_CTW8","Sequence Description", "Enzyme Name", "InterPro Name","Annotation GO ID","Annotation GO Term")
down %<>% mutate(Delta=VSTmedian_CTW8 - VSTmedian_CTW4)
down <- down[,c(1,4:8,2,3,9)]

write_tsv(down,file =  here("../FT1 CRISPR/T89.DEdn.annotated.tsv"))




## How many genes associated to Go category for all, up and down DE genes at CTW8 in FT1 vs T89
#all
all_term_to_genes <- as.data.frame(table(CTW8.T89.FT1_all_enriched_term_to_genes$GOID))
colnames(all_term_to_genes) <- c("GOID", "Number_of_genes")
enrall <- read_tsv(here("../FT1 CRISPR/CTW8.T89.FT1.all-genes_GO-enrichment.tsv"),
                   show_col_types = FALSE) %>% rename(GOID=id)
all_term_to_genes <-  merge(enrall[,c(2,5)], all_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#up
up_term_to_genes <- as.data.frame(table(CTW8.T89.FT1_up_enriched_term_to_genes$GOID))
colnames(up_term_to_genes) <- c("GOID", "Number_of_genes")
enr <- read_tsv(here("../FT1 CRISPR/CTW8.T89.FT1.up-genes_GO-enrichment.tsv"),
                show_col_types = FALSE) %>% rename(GOID=id)
up_term_to_genes <-  merge(enr[,c(2,5)], up_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#down
dn_term_to_genes <- as.data.frame(table(CTW8.T89.FT1_dn_enriched_term_to_genes$GOID))
colnames(dn_term_to_genes) <- c("GOID", "Number_of_genes")
enrdn <- read_tsv(here("../FT1 CRISPR/CTW8.T89.FT1.dn-genes_GO-enrichment.tsv"),
                  show_col_types = FALSE) %>% rename(GOID=id)
dn_term_to_genes <-  merge(enrdn[,c(2,5)], dn_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

## How many genes associated to Go category for all, up and down DE genes at CTW8 vs CTW4 in FT1
 # all
all_term_to_genes <- as.data.frame(table(FT1.CTW8.CTW4.all_enriched_term_to_genes$GOID))
colnames(all_term_to_genes) <- c("GOID", "Number_of_genes")
enrall <- read_tsv(here("../FT1 CRISPR/FT1.CTW8.CTW4.all-genes_GO-enrichment.tsv"),
                   show_col_types = FALSE) %>% rename(GOID=id)
all_term_to_genes <-  merge(enrall[,c(2,5)], all_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#up
up_term_to_genes <- as.data.frame(table(FT1.CTW8.CTW4.up_enriched_term_to_genes$GOID))
colnames(up_term_to_genes) <- c("GOID", "Number_of_genes")
enr <- read_tsv(here("../FT1 CRISPR/FT1.CTW8.CTW4.up-genes_GO-enrichment.tsv"),
                show_col_types = FALSE) %>% rename(GOID=id)
up_term_to_genes <-  merge(enr[,c(2,5)], up_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#down
dn_term_to_genes <- as.data.frame(table(FT1.CTW8.CTW4.dn_enriched_term_to_genes$GOID))
colnames(dn_term_to_genes) <- c("GOID", "Number_of_genes")
enrdn <- read_tsv(here("../FT1 CRISPR/FT1.CTW8.CTW4.dn-genes_GO-enrichment.tsv"),
                  show_col_types = FALSE) %>% rename(GOID=id)
dn_term_to_genes <-  merge(enrdn[,c(2,5)], dn_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

## How many genes associated to Go category for all, up and down DE genes at CTW8 vs CTW4 in T89
#all
all_term_to_genes <- as.data.frame(table(T89.CTW8.CTW4.all_enriched_term_to_genes$GOID))
colnames(all_term_to_genes) <- c("GOID", "Number_of_genes")
enrall <- read_tsv(here("../FT1 CRISPR/T89.CTW8.CTW4.all-genes_GO-enrichment.tsv"),
                   show_col_types = FALSE) %>% rename(GOID=id)
all_term_to_genes <-  merge(enrall[,c(2,5)], all_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#up
up_term_to_genes <- as.data.frame(table(T89.CTW8.CTW4.up_enriched_term_to_genes$GOID))
colnames(up_term_to_genes) <- c("GOID", "Number_of_genes")
enr <- read_tsv(here("../FT1 CRISPR/T89.CTW8.CTW4.up-genes_GO-enrichment.tsv"),
                show_col_types = FALSE) %>% rename(GOID=id)
up_term_to_genes <-  merge(enr[,c(2,5)], up_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#down
dn_term_to_genes <- as.data.frame(table(T89.CTW8.CTW4.dn_enriched_term_to_genes$GOID))
colnames(dn_term_to_genes) <- c("GOID", "Number_of_genes")
enrdn <- read_tsv(here("../FT1 CRISPR/T89.CTW8.CTW4.dn-genes_GO-enrichment.tsv"),
                  show_col_types = FALSE) %>% rename(GOID=id)
dn_term_to_genes <-  merge(enrdn[,c(2,5)], dn_term_to_genes, by.y= "GOID", all.x= TRUE, all.y= TRUE)

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


