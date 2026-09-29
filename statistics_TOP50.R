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

## ---- Empty result tables ----
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
  
  g_patient <- build_sample_network(cormat, toptable_edges, patient_col, n_top = 50)
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
saveRDS(graph_stats, file.path("new_data", "graph_stats50.rds"))
saveRDS(cd44_stats,  file.path("new_data", "cd44_stats50.rds"))


data11 <- readRDS(file.path("new_data", "graph_stats50.rds"))
data22 <- readRDS(file.path("new_data", "cd44_stats50.rds"))

View(data11)
View(data22)
