#!/usr/bin/env Rscript

library(readxl)

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument))
analysis_dir <- dirname(dirname(script_path))
raw_dir <- file.path(analysis_dir, "data", "raw")
table_dir <- file.path(analysis_dir, "results", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

target_row <- function(data, id_column, gene_id, path) {
  index <- which(data[[id_column]] == gene_id)
  if (length(index) != 1L) {
    stop("Expected one row for ", gene_id, " in ", basename(path), ".")
  }
  data[index, , drop = FALSE]
}

read_target_values <- function(path, gene_id) {
  expression_matrix <- read.delim(path, check.names = FALSE)
  target <- target_row(expression_matrix, 1L, gene_id, path)
  values <- as.numeric(target[1L, -1L])
  names(values) <- colnames(expression_matrix)[-1L]
  if (any(!is.finite(values) | values < 0)) {
    stop("Missing or invalid scaled abundance for ", gene_id, ".")
  }
  values
}

extract_developmental <- function(values, species, strain, gene_id) {
  keep <- grepl("^bio[0-9]+\\.hr[0-9]+$", names(values))
  parsed <- strcapture(
    "^bio([0-9]+)\\.hr([0-9]+)$", names(values)[keep],
    proto = list(replicate_number = integer(), time_hour = integer())
  )
  result <- data.frame(
    species = species,
    strain = strain,
    gene_id = gene_id,
    biological_replicate = paste0("bio", parsed$replicate_number),
    time_hour = parsed$time_hour,
    scaled_abundance = unname(values[keep])
  )
  # the +1 keeps zero-abundance observations on plotting scale
  result$log2_scaled_abundance <- log2(result$scaled_abundance + 1)
  result
}

extract_celltype_pairs <- function(values, species, strain, gene_id) {
  sample_names <- names(values)[
    grepl("^bio[0-9]+\\.(prespore|prestalk)$", names(values))
  ]
  replicate_ids <- unique(sub("\\..*$", "", sample_names))
  prespore <- unname(values[paste0(replicate_ids, ".prespore")])
  prestalk <- unname(values[paste0(replicate_ids, ".prestalk")])
  if (length(replicate_ids) == 0L || anyNA(c(prespore, prestalk)) ||
      any(prestalk <= 0)) {
    stop("Missing cell-type pair or zero prestalk denominator for ", species, ".")
  }
  data.frame(
    species = species,
    strain = strain,
    gene_id = gene_id,
    biological_replicate = replicate_ids,
    prestalk_scaled_abundance = prestalk,
    prespore_scaled_abundance = prespore,
    prespore_to_prestalk_ratio = prespore / prestalk
  )
}

read_published_p_value <- function(path, gene_id, id_column) {
  data <- as.data.frame(read_excel(
    path, skip = 18, col_types = c("text", rep("guess", 10))
  ))
  target <- target_row(data, id_column, gene_id, path)

  target$ebays.pval
}

specifications <- list(
  list(
    species = "D. discoideum",
    strain_developmental = "AX4",
    strain_celltype = "NC4",
    gene_id = "DDB_G0278731",
    id_column = "ddb_g",
    scaled_file = "discoideum_normalized_data_set.txt",
    statistics_file = "Supplementary Table 5.xls"
  ),
  list(
    species = "D. purpureum",
    strain_developmental = "DpAX1",
    strain_celltype = "DpAX1",
    gene_id = "87341",
    id_column = "jgi_id",
    scaled_file = "purpureum_normalized_data_set.txt",
    statistics_file = "Supplementary Table 6.xls"
  )
)

developmental_results <- list()
celltype_results <- list()
for (i in seq_along(specifications)) {
  s <- specifications[[i]]
  values <- read_target_values(file.path(raw_dir, s$scaled_file), s$gene_id)
  developmental_results[[i]] <- extract_developmental(
    values, s$species, s$strain_developmental, s$gene_id
  )
  pairs <- extract_celltype_pairs(values, s$species, s$strain_celltype, s$gene_id)
  pairs$whole_transcript_ebayes_p_value <- read_published_p_value(
    file.path(raw_dir, s$statistics_file), s$gene_id, s$id_column
  )
  celltype_results[[i]] <- pairs
}

developmental <- do.call(rbind, developmental_results)
prespore <- do.call(rbind, celltype_results)

# save values behind the panels as reusable analysis outputs
write_tsv(developmental, file.path(table_dir, "parikh_pgtd_developmental_expression.tsv"))
write_tsv(prespore, file.path(table_dir, "parikh_pgtd_prespore_enrichment.tsv"))

# wrapper launches two R sessions so record each stage so readxl is included
writeLines(
  c("Extraction: scripts/01_extract_parikh_pgtd.R", capture.output(sessionInfo())),
  file.path(analysis_dir, "results", "session_info.txt")
)
message("Extracted the Parikh developmental and cell-type pgtD results")
