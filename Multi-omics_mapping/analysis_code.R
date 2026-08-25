# ##############################################################################
#
##  Data analysis script
### For gutcin project
#
# ##############################################################################

# general ##############################################################################

# Packages
# Packages
library(dplyr)
library(ggplot2)
library(DESeq2)
library(tibble)

#work diretory (your diretory)
setwd("./Multi-omics_mapping/")
set.seed(123)

# diversity ##############################################################################
#read data
counts_data <- as.matrix(read.table('00.output/counts.txt', sep='\t', header=TRUE, 
                                    stringsAsFactors = FALSE, 
                                    check.names = FALSE, quote='',row.names = 1))
counts_data1 <- counts_data[,-1]
metadata <- read.csv('00.data/metadata.csv', sep=',')

# Ensure sample names match
common_samples <- intersect(colnames(counts_data), metadata$Samples)
counts_data <- counts_data[, common_samples]
metadata <- metadata[metadata$Samples %in% common_samples, ]

# 1. Calculate the number of genes (expressed genes) in each sample.
gene_counts <- data.frame(
  Samples = colnames(counts_data),
  ExpressedGenes = colSums(counts_data > 0)
)

# 2. Calculate the diversity (Shannon diversity index) for each sample.
shannon_diversity <- function(x) {
  x <- x[x > 0]
  p <- x / sum(x)
  -sum(p * log(p))
}

diversity_data <- data.frame(
  Samples = colnames(counts_data),
  ShannonDiversity = apply(counts_data, 2, shannon_diversity)
)

# Combine
result_data <- metadata %>%
  left_join(gene_counts, by = "Samples") %>%
  left_join(diversity_data, by = "Samples")

write.csv(result_data, "./00.output/diversity.csv", row.names = FALSE)

# deseq2 ##############################################################################
#Data processing and confirmation
expr_matrix <- as.matrix(counts_data1)+1
stopifnot(all(metadata$Samples %in% colnames(expr_matrix)))
studies <- metadata %>% pull(Datasets) %>% unique

# Function: Perform DESeq2 analysis on a single dataset
run_deseq2_per_dataset <- function(expr_matrix, metadata, dataset_name) {
  # Filter samples in the current dataset
  dataset_samples <- metadata$Samples[metadata$Datasets == dataset_name]
  dataset_expr <- expr_matrix[, colnames(expr_matrix) %in% dataset_samples]
  dataset_meta <- metadata[metadata$Datasets == dataset_name, ]
  
  # Ensure sample order is consistent
  dataset_meta <- dataset_meta[match(colnames(dataset_expr), dataset_meta$Samples), ]
  
  # Create a DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = dataset_expr,
    colData = dataset_meta,
    design = ~ Group
  )
  nrow(dds)
  # Run DESeq2 analysis
  dds$Group_new <- factor(ifelse(dds$Group %in% c("CRC", "UC"), "Disease", "Healthy"))
  design(dds) <- ~ Group_new
  
  dds <- DESeq(dds)
  
  # Extraction Results
  res <- results(dds, contrast = c("Group_new", "Healthy", "Disease"))
  
  # Convert to a data frame and add necessary information
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      dataset = dataset_name,
      n_total = ncol(dataset_expr),
      n_UC = sum(dataset_meta$Group == "Disease"),
      n_Healthy = sum(dataset_meta$Group == "Healthy")
    )
  
  return(res_df)
}

# Main function
run_cross_dataset_analysis <- function(expr_matrix, metadata) {
  
  datasets <- unique(metadata$Datasets)
  cat("Find", length(datasets), "Datasets:", paste(datasets, collapse = ", "), "\n")
  
  all_results <- list()
  
  for (dataset in datasets) {
    cat("In analysis:", dataset, "\n")
    
    tryCatch({
      dataset_result <- run_deseq2_per_dataset(expr_matrix, metadata, dataset)
      all_results[[dataset]] <- dataset_result
      cat("  Finished! Obtained", nrow(dataset_result), "gene\n")
    }, error = function(e) {
      cat("  Error analyzing", dataset, ":", e$message, "\n")
    })
  }
  
  # Combine result
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  
  cat("Analysis complete!\n")
  
  return(combined_results)
}

# Run analysis
results <- run_cross_dataset_analysis(expr_matrix, metadata)

# Save result
write.csv(results, "./00.output/individual_dataset_results.csv", row.names = FALSE)
