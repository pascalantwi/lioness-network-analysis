

dir.create(file.path("output", "figures"), showWarnings = FALSE, recursive = TRUE)


##Read in preprocessed data 
exp     <- readRDS(file.path("analysis", "exp.rds"))
targets <- readRDS(file.path("analysis", "targets.rds"))


##Pathway database 
msigdb_data   <- getMsigdb(org = 'hs', id = 'SYM', version = '2023.1')
reactome_sets <- subsetCollection(msigdb_data, subcollection = 'CP:REACTOME')



#Pull the genes for one Reactome pathway
get_pathway_genes <- function(pathway_name, reactome_sets, exp) {
  stopifnot(pathway_name %in% names(reactome_sets))
  pathway_genes <- geneIds(reactome_sets[[pathway_name]])
  intersect(pathway_genes, rownames(exp))
}

# Identify edges whose |correlation difference| exceeds a threshold
select_edges <- function(netdiff, threshold) {
  melted <- melt(upper.tri(netdiff))
  melted <- melted[which(melted$value), ]
  values <- netdiff[which(upper.tri(netdiff))]
  melted <- cbind(melted[, 1:2], values)
  gene_names <- row.names(netdiff)
  melted[,1] <- gene_names[melted[,1]]
  melted[,2] <- gene_names[melted[,2]]
  row.names(melted) <- paste(melted[,1], melted[,2], sep = "_")
  row.names(melted[which(abs(melted[,3]) > threshold), ])
}

#The full LIONESS + LIMMA + network plot pipeline on one gene set
#The LIMMA results table, and writes one PNG figure to output_file
run_lioness_analysis <- function(exp, targets, genes, output_file,
                                 edge_threshold = 0.5, top_n_edges = 50) {
  
  dat <- as.matrix(exp[genes, ])
  
  groupyes <- which(targets$mets == "yes")
  groupno  <- which(targets$mets == "no")
  netyes   <- cor(t(dat[, groupyes]))
  netno    <- cor(t(dat[, groupno]))
  netdiff  <- netyes - netno
  
  tosel <- select_edges(netdiff, threshold = edge_threshold)
  
  cormat <- lioness(dat, netFun)
  row.names(cormat) <- paste(cormat[,1], cormat[,2], sep = "_")
  corsub <- cormat[which(row.names(cormat) %in% tosel), 3:ncol(cormat)]
  corsub <- as.matrix(corsub)
  
  group <- factor(targets$mets)
  design <- model.matrix(~0 + group)
  cont.matrix <- makeContrasts(yesvsno = groupyes - groupno, levels = design)
  fit    <- lmFit(corsub, design)
  fit2   <- contrasts.fit(fit, cont.matrix)
  fit2e  <- eBayes(fit2)
  toptable <- topTable(fit2e, number = nrow(corsub), adjust = "fdr")
  
  n_top <- min(top_n_edges, nrow(toptable))
  toptable_edges <- t(matrix(unlist(strsplit(row.names(toptable)[1:n_top], "_")), 2))
  z <- cbind(toptable_edges, toptable$logFC[1:n_top])
  g <- graph.data.frame(z, directed = FALSE)
  E(g)$weight <- as.numeric(z[,3])
  E(g)$color[E(g)$weight < 0] <- "blue"
  E(g)$color[E(g)$weight > 0] <- "red"
  E(g)$weight <- 1
  
  topgeneslist <- unique(c(toptable_edges))
  fit_de   <- lmFit(exp, design)
  fit2_de  <- contrasts.fit(fit_de, cont.matrix)
  fit2e_de <- eBayes(fit2_de)
  topDE    <- topTable(fit2e_de, number = nrow(exp), adjust = "fdr")
  topDE    <- topDE[which(row.names(topDE) %in% topgeneslist), ]
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
  bincol   <- sapply(V(g)$weight, function(w) min(which(breaks4 > w)))
  V(g)$color <- mypalette4[bincol]
  
  png(output_file, width = 1200, height = 1200, res = 150)
  par(mar = c(0,0,0,0))
  plot(g, vertex.label.cex = 0.7, vertex.size = 10, vertex.label.color = "black",
       vertex.label.font = 3, edge.width = 10*(abs(as.numeric(z[,3])) - 0.7),
       vertex.color = V(g)$color)
  dev.off()
  
  toptable
}



##genes for a pathway

get_pathway_genes <- function(pathway_name, reactome_sets, exp) {
  if (pathway_name == "CUSTOM_IMMUNE_SYSTEM_UNION") {
    sub_names <- grep("IMMUNE_SYSTEM", names(reactome_sets), value = TRUE)
    pathway_genes <- unique(unlist(lapply(sub_names, function(nm) geneIds(reactome_sets[[nm]]))))
  } else {
    stopifnot(pathway_name %in% names(reactome_sets))
    pathway_genes <- geneIds(reactome_sets[[pathway_name]])
  }
  intersect(pathway_genes, rownames(exp))
}


##PATHWAYS 

pathways <- data.frame(
  reactome_id  = c("R-HSA-1280218", "R-HSA-1280215", "R-HSA-913531", "R-HSA-168256"),
  pathway_name = c("REACTOME_ADAPTIVE_IMMUNE_SYSTEM",
                   "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
                   "REACTOME_INTERFERON_SIGNALING",
                   "CUSTOM_IMMUNE_SYSTEM_UNION"),
  stringsAsFactors = FALSE
)


pathways


## all  pathways loop 

pathway_results <- list()

for (i in seq_len(nrow(pathways))) {
  
  reactome_id  <- pathways$reactome_id[i]
  pathway_name <- pathways$pathway_name[i]
  
  pathway_genes <- get_pathway_genes(pathway_name, reactome_sets, exp)
  cat(reactome_id, "-", pathway_name, "-", length(pathway_genes), "genes present in exp\n")
  
  output_file <- file.path("output", "figures",
                           paste0(reactome_id, "_", pathway_name, "_network.png"))
  
  pathway_results[[reactome_id]] <- run_lioness_analysis(
    exp, targets, pathway_genes, output_file,
    edge_threshold = 0.5, top_n_edges = 50
  )
}


##output 
list.files(file.path("output", "figures"))
lapply(pathway_results, head)


