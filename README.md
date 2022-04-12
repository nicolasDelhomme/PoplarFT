# FLOWERING LOCUS T Paralogs Control the Annual Growth Cycle in Populus Trees

**Domenique André1, Keh Chien Lee1, Daniela Goretti1¶, Alice Marcon1, Bo Zhang1, Nicolas Delhomme1, Markus Schmid2 and Ove Nilsson1+**

1 Umeå Plant Science Centre, Department of Forest Genetics and Plant Physiology, Swedish University of Agricultural Sciences, 901 83 Umeå, Sweden

2 Umeå Plant Science Centre, Department of Plant Physiology, Umeå University, 907 36 Umeå, Sweden
¶Present address: Umeå Plant Science Centre, Department of Plant Physiology, Umeå University, 907 36 Umeå, Sweden

+ Corresponding author: Ove Nilsson (ove.nilsson@slu.se)

## Content
This repository contains the code used to perform the preprocessing and the analysis of the data. 

1. All the code for the preprocessing is available in the UPSCb-common submodule pipeline repository. 

2. The results of the preprocessing has been summarised using MultiQC and is available in the [doc](https://github.com/nicolasDelhomme/PoplarFT/tree/main/doc) directory as an html report.

3. The code for the data analysis, both the exploratory data analysis (EDA) and differential expression (DE) analysis are available in the [src/R](https://github.com/nicolasDelhomme/PoplarFT/tree/main/src/R) directory. The EDA is in the BiologicalQA.R file, while the DE is in the DifferentialExpression.R file.

4. The metadata information for these analysis is available in the [doc](https://github.com/nicolasDelhomme/PoplarFT/tree/main/doc) directory, in the csv file.

## Summary
In temperate and boreal regions, perennials adapt their annual growth cycle to the change of seasons. These adaptations ensure survival in harsh environmental conditions, allowing growth at different latitudes and altitudes, and are therefore tightly regulated. Populus tree species cease growth and form terminal buds in autumn when photoperiod falls below a certain threshold [1]. This is followed by establishment of dormancy and cold hardiness over the winter. At the center of the photoperiodic pathway in Populus is the gene FLOWERING LOCUS T2 (FT2), which is expressed during summer and harbors significant SNPs in its locus associated with timing of bud set [1-4]. The paralogous gene FT1, on the other hand, is hyper-induced in chilling buds during winter [3, 5]. Even though its function is so far unknown, it has been suggested to be involved in the regulation of flowering and the release of winter dormancy [3, 5]. 
In this study we employ CRISPR/Cas9-mediated gene editing to individually study the function of the FT-like genes in Populus trees. We show that while FT2 is required for vegetative growth during spring and summer and regulates the entry into dormancy, expression of FT1 is absolutely required for bud flush in spring. Gene expression profiling suggests that this function of FT1 is linked to the release of winter dormancy rather than to the regulation of bud flush per se.
These data show how FT duplication and sub-functionalization have allowed Populus trees to regulate two completely different and major developmental control points during the yearly growth cycle.

## Keywords
Populus, FLOWERING LOCUS T, Paralogs, Annual growth cycle, Dormancy, Bud flush
