### METRIC

build_full_sample_network <- function(cormat, patient_col) {
  gene1 <- as.character(cormat[, 1])
  gene2 <- as.character(cormat[, 2])
  
  ## keep only ONE row per unique gene pair
  keep <- gene1 < gene2
  
  edge_list <- data.frame(
    gene1  = gene1[keep],
    gene2  = gene2[keep],
    weight = as.numeric(cormat[keep, patient_col])
  )
  
  g_patient <- graph_from_data_frame(edge_list, directed = FALSE)
  E(g_patient)$weight <- edge_list$weight
  g_patient
}

## graph-level stats for one patient's network
get_graph_stats <- function(g_patient, patient_col) {
  data.frame(
    patient            = patient_col,
    n_nodes            = vcount(g_patient),
    n_edges            = ecount(g_patient),
    diameter           = igraph::diameter(g_patient, weights = NA),
    transitivity       = igraph::transitivity(g_patient, type = "global"),
    total_edge_weight  = sum(E(g_patient)$abs_weight)
  )
}

## CD44 node-level
get_cd44_stats <- function(g_patient, patient_col) {
  if (!("CD44" %in% V(g_patient)$name)) {
    return(NULL)
  }
  
  data.frame(
    patient      = patient_col,
    degree       = igraph::strength(g_patient, v = "CD44", weights = E(g_patient)$abs_weight),
    betweenness  = igraph::betweenness(g_patient, v = "CD44", weights = E(g_patient)$abs_weight)
  )
}

## patient ID -> mets status ("yes"/"no")
mets_lookup <- setNames(targets$mets, targets$sample)

## Empty result tables
graph_stats <- data.frame(
  patient = character(), n_nodes = integer(), n_edges = integer(),
  diameter = numeric(), transitivity = numeric(),
  total_edge_weight = numeric(), mets = character(),
  stringsAsFactors = FALSE
)

cd44_stats <- data.frame(
  patient = character(),
  degree = numeric(), betweenness = numeric(), mets = character(),
  stringsAsFactors = FALSE
)

## Loop over EVERY patient
for (patient_col in patient_columns) {
  
  cat("\n===== Patient:", patient_col, "=====\n")
  
  g_patient <- build_full_sample_network(cormat, patient_col)
  E(g_patient)$abs_weight <- abs(E(g_patient)$weight)
  
  graph_row <- get_graph_stats(g_patient, patient_col)
  graph_row$mets <- mets_lookup[[patient_col]]
  cat("Graph-level stats for", patient_col, "(mets =", graph_row$mets, "):\n")
  print(graph_row)
  graph_stats <- rbind(graph_stats, graph_row)
  
  cd44_row <- get_cd44_stats(g_patient, patient_col)
  if (!is.null(cd44_row)) {
    cd44_row$mets <- mets_lookup[[patient_col]]
    cat("CD44 stats for", patient_col, "(mets =", cd44_row$mets, "):\n")
    print(cd44_row)
    cd44_stats <- rbind(cd44_stats, cd44_row)
  } else {
    cat("CD44 not present in", patient_col, "'s network — skipped\n")
  }
}

## Save results
dir.create("new_data", showWarnings = FALSE, recursive = TRUE)
saveRDS(graph_stats, file.path("new_data", "graph_stats_full.rds"))
saveRDS(cd44_stats,  file.path("new_data", "cd44_stats_full.rds"))

data1 <- readRDS(file.path("new_data", "graph_stats_full.rds"))
data2 <- readRDS(file.path("new_data", "cd44_stats_full.rds"))

View(data1)
View(data2)
