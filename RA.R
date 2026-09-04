
# Data from Geo
load("C:/Users/HP/Desktop/OSdata.RData")


#Packages 
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "lumi"))
BiocManager::install("illuminaHumanv2.db")

# Libraries 
library(AnnotationDbi)
library(illuminaHumanv2.db)
library(GEOquery)

#Data 
data_raw <- getGEO(filename = "C:/Users/HP/OneDrive/Desktop/RA/GSE42352_series/GSE42352_series_matrix.txt")
class(data_raw)



#exp_all and sample_info
exp_all     <- exprs(data_raw)
sample_info <- pData(data_raw)

dim(exp_all)          #genes x samples
dim(sample_info)       #samples x metadata fields
rownames(exp_all)[1:5] #probe/nuID names
colnames(exp_all)[1:5] #sample GSM accessions
head(sample_info$title)


#samples and the outcome variable 
colnames(sample_info)               
head(sample_info$title, 20)         

### Filtering 
# 84 patient biopsies 
is_biopsy   <- grepl("biopsy", sample_info$title)
sum(is_biopsy)

biopsy_exp  <- exp_all[, is_biopsy]
biopsy_meta <- sample_info[is_biopsy, ]
dim(biopsy_exp)


#lplatform annotation table 
gpl            <- getGEO("GPL10295")
platform_table <- Table(gpl)

colnames(platform_table)
head(platform_table)


#map probes to gene symbols 
id_col     <- colnames(platform_table)[1]
symbol_col <- colnames(platform_table)[grepl("symbol", colnames(platform_table), ignore.case = TRUE)][1]
id_col
symbol_col

match_pos  <- match(rownames(biopsy_exp), platform_table[[id_col]])
gene_names <- platform_table[[symbol_col]][match_pos]

sum(is.na(gene_names))   #probes with no match, get dropped below


biopsy_data      <- as.data.frame(biopsy_exp)
biopsy_data$gene <- gene_names
biopsy_data      <- biopsy_data[!is.na(biopsy_data$gene) & biopsy_data$gene != "", ]


###: Collapsing duplicate genes

#collapse duplicate gene symbols
sample_cols        <- setdiff(colnames(biopsy_data), "gene")
biopsy_data$spread <- apply(biopsy_data[, sample_cols], 1, var)
biopsy_data        <- biopsy_data[order(biopsy_data$gene, -biopsy_data$spread), ]
biopsy_data        <- biopsy_data[!duplicated(biopsy_data$gene), ]
rownames(biopsy_data) <- biopsy_data$gene

full_exp <- as.matrix(biopsy_data[, sample_cols])
dim(full_exp)   

## Building final


# build exp and targets, subset to patients with known outcome
mets_status   <- biopsy_meta[["metastasis within 5yrs:ch1"]]
table(mets_status)

has_mets_info <- !is.na(mets_status) & mets_status != ""
sum(has_mets_info)   

exp     <- full_exp[, has_mets_info]
targets <- data.frame(sample = colnames(exp), mets = mets_status[has_mets_info])
rownames(targets) <- targets$sample

dim(exp)
dim(targets)
table(targets$mets)



#process data 
#variables 
str(exp)             
str(targets)          

head(exp)             
head(targets)       

dim(exp)
dim(targets)


rownames(exp)[1:10]    
colnames(exp)[1:10]    
colnames(targets) 


#mean and sd 
mean(exp)                     
sd(exp)                       

#mean and sd per patient group
tapply(exp, rep(targets$mets, each = nrow(exp)), mean)
tapply(exp, rep(targets$mets, each = nrow(exp)), sd)

# percentages
table(targets$mets)                  
round(prop.table(table(targets$mets)) * 100, 1)   

# ACCROSS
range(exp)          

#range per patient group
range(exp[, targets$sample[targets$mets == "yes"]])
range(exp[, targets$sample[targets$mets == "no"]])


#### APPLICATION OF THE LIONESS

install.packages("devtools")
devtools::install_github("kuijjerlab/lionessR")
library(lionessR)
library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)

#subset to the 500 most variable genes 
nsel <- 500
cvar <- apply(as.array(as.matrix(exp)), 1, sd)
dat  <- cbind(cvar, exp)
dat  <- dat[order(dat[,1], decreasing = TRUE), ]
dat  <- dat[1:nsel, -1]
dat  <- as.matrix(dat)

#two condition-specific networks, and their difference 
groupyes <- which(targets$mets == "yes")
groupno  <- which(targets$mets == "no")
netyes   <- cor(t(dat[, groupyes]))
netno    <- cor(t(dat[, groupno]))
netdiff  <- netyes - netno

# convert to an edge list, keep edges with |difference| > 0.5 
cormat2 <- rep(1:nsel, each = nsel)
cormat1 <- rep(1:nsel, nsel)
el      <- cbind(cormat1, cormat2, c(netdiff))
melted  <- melt(upper.tri(netdiff))
melted  <- melted[which(melted$value), ]
values  <- netdiff[which(upper.tri(netdiff))]
melted  <- cbind(melted[, 1:2], values)
genes   <- row.names(netdiff)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep = "_")
tosub <- melted
tosel <- row.names(tosub[which(abs(tosub[,3]) > 0.5), ])

# run LIONESS on the full 500-gene set, subset to those edges 
cormat <- lioness(dat, netFun)
row.names(cormat) <- paste(cormat[,1], cormat[,2], sep = "_")
corsub <- cormat[which(row.names(cormat) %in% tosel), 3:ncol(cormat)]
corsub <- as.matrix(corsub)

# LIMMA: differential co-expression between groups 
group <- factor(targets$mets)
design <- model.matrix(~0 + group)
cont.matrix <- makeContrasts(yesvsno = (groupyes - groupno), levels = design)
fit    <- lmFit(corsub, design)
fit2   <- contrasts.fit(fit, cont.matrix)
fit2e  <- eBayes(fit2)
toptable <- topTable(fit2e, number = nrow(corsub), adjust = "fdr")

head(toptable)   


#50 differential edges as a network graph 
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(toptable), "_"))), 2))
z <- cbind(toptable_edges[1:50, ], toptable$logFC[1:50])
g <- graph.data.frame(z, directed = FALSE)
E(g)$weight <- as.numeric(z[,3])
E(g)$color[E(g)$weight < 0] <- "blue"
E(g)$color[E(g)$weight > 0] <- "red"
E(g)$weight <- 1

#differential expression, to color nodes 
topgeneslist <- unique(c(toptable_edges[1:50, ]))
fit   <- lmFit(exp, design)
fit2  <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
topDE <- topTable(fit2e, number = nrow(exp), adjust = "fdr")
topDE <- topDE[which(row.names(topDE) %in% topgeneslist), ]
topgenesDE <- cbind(row.names(topDE), topDE$t)

nodeorder <- cbind(V(g)$name, 1:length(V(g)))
nodes <- merge(nodeorder, topgenesDE, by.x = 1, by.y = 1)
nodes <- nodes[order(as.numeric(as.character(nodes[,2]))), ]
nodes[,3] <- as.numeric(as.character(nodes[,3]))
nodes <- nodes[, -2]
V(g)$weight <- nodes[,2]

mypalette4 <- colorRampPalette(c("blue","white","white","red"), space = "Lab")(256)
breaks2a <- seq(min(V(g)$weight), 0, length.out = 128)
breaks2b <- seq(0.00001, max(V(g)$weight) + 0.1, length.out = 128)
breaks4  <- c(breaks2a, breaks2b)

bincol <- rep(NA, length(V(g)))
for (i in 1:length(V(g))) {
  bincol[i] <- min(which(breaks4 > V(g)$weight[i]))
}
bincol <- mypalette4[bincol]
V(g)$color <- bincol


# plot
png("C:/Users/HP/OneDrive/Desktop/RA/lioness_network.png",
    width = 1200, height = 1200, res = 150)

par(mar = c(0,0,0,0))
plot(g, vertex.label.cex = 0.7, vertex.size = 10, vertex.label.color = "black",
     vertex.label.font = 3, edge.width = 10*(abs(as.numeric(z[,3])) - 0.7),
     vertex.color = V(g)$color)

dev.off()



## Different subset

nsel        <- 1000   
edge_thresh <- 0.5   
top_n_edges <- 100    

#subset to top nsel most variable genes 
cvar <- apply(as.array(as.matrix(exp)), 1, sd)
dat  <- cbind(cvar, exp)
dat  <- dat[order(dat[,1], decreasing = TRUE), ]
dat  <- dat[1:nsel, -1]
dat  <- as.matrix(dat)

#two condition-specific networks, their difference 
groupyes <- which(targets$mets == "yes")
groupno  <- which(targets$mets == "no")
netyes   <- cor(t(dat[, groupyes]))
netno    <- cor(t(dat[, groupno]))
netdiff  <- netyes - netno

#edge list, keep edges above edge_thresh 
melted  <- melt(upper.tri(netdiff))
melted  <- melted[which(melted$value), ]
values  <- netdiff[which(upper.tri(netdiff))]
melted  <- cbind(melted[, 1:2], values)
genes   <- row.names(netdiff)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep = "_")
tosel <- row.names(melted[which(abs(melted[,3]) > edge_thresh), ])
length(tosel)   

# LIONESS on the nsel-gene set
cormat <- lioness(dat, netFun)
row.names(cormat) <- paste(cormat[,1], cormat[,2], sep = "_")
corsub <- cormat[which(row.names(cormat) %in% tosel), 3:ncol(cormat)]
corsub <- as.matrix(corsub)

# LIMMA differential co-expression 
group <- factor(targets$mets)
design <- model.matrix(~0 + group)
cont.matrix <- makeContrasts(yesvsno = (groupyes - groupno), levels = design)
fit    <- lmFit(corsub, design)
fit2   <- contrasts.fit(fit, cont.matrix)
fit2e  <- eBayes(fit2)
toptable <- topTable(fit2e, number = nrow(corsub), adjust = "fdr")
head(toptable)
#take top_n_edges, build network graph 
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(toptable), "_"))), 2))
z <- cbind(toptable_edges[1:top_n_edges, ], toptable$logFC[1:top_n_edges])
g <- graph.data.frame(z, directed = FALSE)
E(g)$weight <- as.numeric(z[,3])
E(g)$color[E(g)$weight < 0] <- "blue"
E(g)$color[E(g)$weight > 0] <- "red"
E(g)$weight <- 1

#node coloring by differential expression 
topgeneslist <- unique(c(toptable_edges[1:top_n_edges, ]))
fit   <- lmFit(exp, design)
fit2  <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
topDE <- topTable(fit2e, number = nrow(exp), adjust = "fdr")
topDE <- topDE[which(row.names(topDE) %in% topgeneslist), ]
topgenesDE <- cbind(row.names(topDE), topDE$t)

nodeorder <- cbind(V(g)$name, 1:length(V(g)))
nodes <- merge(nodeorder, topgenesDE, by.x = 1, by.y = 1)
nodes <- nodes[order(as.numeric(as.character(nodes[,2]))), ]
nodes[,3] <- as.numeric(as.character(nodes[,3]))
nodes <- nodes[, -2]
V(g)$weight <- nodes[,2]

mypalette4 <- colorRampPalette(c("blue","white","white","red"), space = "Lab")(256)
breaks2a <- seq(min(V(g)$weight), 0, length.out = 128)
breaks2b <- seq(0.00001, max(V(g)$weight) + 0.1, length.out = 128)
breaks4  <- c(breaks2a, breaks2b)
bincol <- sapply(V(g)$weight, function(w) min(which(breaks4 > w)))
V(g)$color <- mypalette4[bincol]

#plot
fname <- paste0("C:/Users/HP/OneDrive/Desktop/RA/lioness_n", nsel,
                "_t", edge_thresh, "_top", top_n_edges, ".png")
png(fname, width = 1200, height = 1200, res = 150)
par(mar = c(0,0,0,0))
plot(g, vertex.label.cex = 0.7, vertex.size = 10, vertex.label.color = "black",
     vertex.label.font = 3, edge.width = 10*(abs(as.numeric(z[,3])) - 0.7),
     vertex.color = V(g)$color)
dev.off()



# Pearson correlation matrices  
dim(netyes)      
dim(netno)

netyes[1:5, 1:5]     #correlation values "yes" group
netno[1:5, 1:5]       #actual correlation values "no" group




install.packages("pheatmap")
library(pheatmap)

#---- same similarity matrix as before ----
network_similarity <- function(cormat) {
  edge_weights <- as.matrix(cormat[, 3:ncol(cormat)])
  cor(edge_weights, method = "spearman")
}

sim_pearson <- network_similarity(cormat)

#Z-score, matching the reference figure's color scale label 
sim_z <- (sim_pearson - mean(sim_pearson)) / sd(sim_pearson)

#group patients together 
group_order <- order(targets$mets)
ann <- data.frame(mets = targets$mets)
rownames(ann) <- targets$sample

#rainbow-style diverging palette
jet_colors <- colorRampPalette(c("darkblue","blue","cyan","green",
                                 "yellow","orange","red"))(100)

pheatmap(sim_z[group_order, group_order],
         cluster_rows = FALSE, cluster_cols = FALSE,
         annotation_col = ann, annotation_row = ann,
         color = jet_colors,
         breaks = seq(-3, 3, length.out = 101),  
         border_color = NA,
         main = "Pearson - Similarity Between Networks",
         filename = "C:/Users/HP/OneDrive/Desktop/RA/pearson_similarity_zscore.png",
         width = 14, height = 13)






#### GITHUB DATA 
load("C:/Users/HP/Downloads/OSdata.RData")

library(reshape2)
library(lionessR)
library(limma)
library(igraph)

rm(toptable)  

nsel <- 500
cvar <- apply(as.array(as.matrix(exp)), 1, sd)
dat  <- cbind(cvar, exp)
dat  <- dat[order(dat[,1], decreasing = TRUE), ]
dat  <- dat[1:nsel, -1]
dat  <- as.matrix(dat)

groupyes <- which(targets$mets == "yes")
groupno  <- which(targets$mets == "no")
netyes   <- cor(t(dat[, groupyes]))
netno    <- cor(t(dat[, groupno]))
netdiff  <- netyes - netno

melted  <- melt(upper.tri(netdiff))
melted  <- melted[which(melted$value), ]
values  <- netdiff[which(upper.tri(netdiff))]
melted  <- cbind(melted[, 1:2], values)
genes   <- row.names(netdiff)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep = "_")
tosel <- row.names(melted[which(abs(melted[,3]) > 0.5), ])

cormat <- lioness(dat, netFun)
row.names(cormat) <- paste(cormat[,1], cormat[,2], sep = "_")
corsub <- cormat[which(row.names(cormat) %in% tosel), 3:ncol(cormat)]
corsub <- as.matrix(corsub)

group <- factor(targets$mets)
design <- model.matrix(~0 + group)
cont.matrix <- makeContrasts(yesvsno = (groupyes - groupno), levels = design)
fit    <- lmFit(corsub, design)
fit2   <- contrasts.fit(fit, cont.matrix)
fit2e  <- eBayes(fit2)
toptable <- topTable(fit2e, number = nrow(corsub), adjust = "fdr")

head(toptable)

# differential edges as a network graph 
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(toptable), "_"))), 2))
z <- cbind(toptable_edges[1:50, ], toptable$logFC[1:50])
g <- graph.data.frame(z, directed = FALSE)
E(g)$weight <- as.numeric(z[,3])
E(g)$color[E(g)$weight < 0] <- "blue"
E(g)$color[E(g)$weight > 0] <- "red"
E(g)$weight <- 1

#differential expression
topgeneslist <- unique(c(toptable_edges[1:50, ]))
fit   <- lmFit(exp, design)
fit2  <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
topDE <- topTable(fit2e, number = nrow(exp), adjust = "fdr")
topDE <- topDE[which(row.names(topDE) %in% topgeneslist), ]
topgenesDE <- cbind(row.names(topDE), topDE$t)

nodeorder <- cbind(V(g)$name, 1:length(V(g)))
nodes <- merge(nodeorder, topgenesDE, by.x = 1, by.y = 1)
nodes <- nodes[order(as.numeric(as.character(nodes[,2]))), ]
nodes[,3] <- as.numeric(as.character(nodes[,3]))
nodes <- nodes[, -2]
V(g)$weight <- nodes[,2]

mypalette4 <- colorRampPalette(c("blue","white","white","red"), space = "Lab")(256)
breaks2a <- seq(min(V(g)$weight), 0, length.out = 128)
breaks2b <- seq(0.00001, max(V(g)$weight) + 0.1, length.out = 128)
breaks4  <- c(breaks2a, breaks2b)

bincol <- rep(NA, length(V(g)))
for (i in 1:length(V(g))) {
  bincol[i] <- min(which(breaks4 > V(g)$weight[i]))
}
bincol <- mypalette4[bincol]
V(g)$color <- bincol

#plot
png("C:/Users/HP/Downloads/OSdata_network.png", width = 1200, height = 1200, res = 150)
par(mar = c(0,0,0,0))
plot(g, vertex.label.cex = 0.7, vertex.size = 10, vertex.label.color = "black",
     vertex.label.font = 3, edge.width = 10*(abs(as.numeric(z[,3])) - 0.7),
     vertex.color = V(g)$color)
dev.off()





