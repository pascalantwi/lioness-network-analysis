#Packages 
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "lumi"))
BiocManager::install("illuminaHumanv2.db")
install.packages("devtools")
install.packages("pheatmap")
devtools::install_github("kuijjerlab/lionessR")
BiocManager::install(c("msigdb", "ExperimentHub", "GSEABase", "clusterProfiler"))


# Libraries 
library(AnnotationDbi)
library(illuminaHumanv2.db)
library(GEOquery)
library(lionessR)
library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)
library(pheatmap)
library(igraph)
library(msigdb)
library(ExperimentHub)
library(GSEABase)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)

