## PREPROCESSING SCRIPT


## Set up folders 
dir.create("data",     showWarnings = FALSE, recursive = TRUE)
dir.create("analysis", showWarnings = FALSE, recursive = TRUE)


## Load data
data_raw_path <- file.path("data", "GSE42352_data_raw.rds")

if (file.exists(data_raw_path)) {
  data_raw <- readRDS(data_raw_path)
} else {
  data_raw <- getGEO("GSE42352", GSEMatrix = TRUE)[[1]]
  saveRDS(data_raw, data_raw_path)
}

class(data_raw)

exp_all     <- exprs(data_raw)
sample_info <- pData(data_raw)

dim(exp_all)
dim(sample_info)


## Filter to the 84 patient biopsies 
is_biopsy   <- grepl("biopsy", sample_info$title)
sum(is_biopsy)

biopsy_exp  <- exp_all[, is_biopsy]
biopsy_meta <- sample_info[is_biopsy, ]
dim(biopsy_exp)


## Platform annotation table
gpl_path <- file.path("data", "GPL10295_platform.rds")

if (file.exists(gpl_path)) {
  platform_table <- readRDS(gpl_path)
} else {
  gpl            <- getGEO("GPL10295")
  platform_table <- Table(gpl)
  saveRDS(platform_table, gpl_path)
}

id_col     <- colnames(platform_table)[1]
symbol_col <- colnames(platform_table)[grepl("symbol", colnames(platform_table), ignore.case = TRUE)][1]


## Map probes to gene symbols 
match_pos  <- match(rownames(biopsy_exp), platform_table[[id_col]])
gene_names <- platform_table[[symbol_col]][match_pos]

biopsy_data      <- as.data.frame(biopsy_exp)
biopsy_data$gene <- gene_names
biopsy_data      <- biopsy_data[!is.na(biopsy_data$gene) & biopsy_data$gene != "", ]


## ---- Collapse duplicate gene symbols, keep highest-variance probe ----
sample_cols        <- setdiff(colnames(biopsy_data), "gene")
biopsy_data$spread <- apply(biopsy_data[, sample_cols], 1, var)
biopsy_data        <- biopsy_data[order(biopsy_data$gene, -biopsy_data$spread), ]
biopsy_data        <- biopsy_data[!duplicated(biopsy_data$gene), ]
rownames(biopsy_data) <- biopsy_data$gene

full_exp <- as.matrix(biopsy_data[, sample_cols])
dim(full_exp)


##  Build final exp/targets
mets_status   <- biopsy_meta[["metastasis within 5yrs:ch1"]]
has_mets_info <- !is.na(mets_status) & mets_status != ""

exp     <- full_exp[, has_mets_info]
targets <- data.frame(sample = colnames(exp), mets = mets_status[has_mets_info])
rownames(targets) <- targets$sample

dim(exp)
dim(targets)
table(targets$mets)


##  Save processed data 
saveRDS(exp,     file.path("analysis", "exp.rds"))
saveRDS(targets, file.path("analysis", "targets.rds"))
