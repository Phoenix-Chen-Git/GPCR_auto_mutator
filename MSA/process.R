library(Biostrings)
library(msa)
library(future.apply)

# Set up parallel plan for 28 CPUs
plan(multisession, workers = 28)

# Function to calculate conservation score (per site) in parallel
conservation_scores <- function(alignment) {
  aln_matrix <- as.matrix(alignment)
  scores <- future_apply(aln_matrix, 2, function(column) {
    freqs <- table(column[column != "-"])
    if (length(freqs) == 0) return(0)
    max(freqs) / sum(freqs)
  })
  return(scores)
}

# Example usage
fasta_file <- "uniprotkb_family_G_protein_coupled_rece_2025_05_25.fasta.gz"
seqs <- readAAStringSet(fasta_file)
seqs_filtered <- seqs[width(seqs) > 250]

# Use ClustalOmega with 32 threads for alignment
alignment <- msa(seqs_filtered, method = "ClustalOmega", 
                 order = "input", param = "--threads 32")

scores <- conservation_scores(alignment)
print(scores)

# Save all alignment conservation scores
write.csv(scores, file = "conservation_scores_all_columns.csv", row.names = FALSE)

# To map a query sequence onto the alignment and get scores for its residues:
query <- "MVNSTHRGMHTSLHLWNRSSYRLHSNASESLGKGYSDGGCYEQLFVSPEVFVTLGVISLLENILVIVAIAKNKNLHSPMYFFICSLAVADMLVSVSNGSETIVITLLNSTDTDAQSFTVNIDNVIDSVICSSLLASICSLLSIAVDRYFTIFYALQYHNIMTVKRVGIIISCIWAACTVSGILFIIYSDSSAVIICLITMFFTMLALMASLYVHMFLMARLHIKRIAVLPGTGAIRQGANMKGAITLTILIGVFVVCWAPFFLHLIFYISCPQNPYCVCFMSHFNLYLILIMCNSIIDPLIYALRSQELRKTFKEIICCYPLGGLCDLSSRY"  # Input sequence as string
input_seq <- AAString(query)

aln_matrix <- as.matrix(alignment)
aln_names <- names(seqs_filtered)
query_idx <- which(sapply(aln_names, function(nm) grepl(query, as.character(seqs_filtered[nm]))))
if (length(query_idx) == 1) {
  query_row <- aln_matrix[query_idx, ]
  query_scores <- scores[query_row != "-"]
  print(query_scores)
  write.csv(query_scores, file = "conservation_scores_query.csv", row.names = FALSE)
} else {
  cat("Query sequence not found in alignment. Manual mapping needed.\n")
}
