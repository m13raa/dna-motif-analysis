library (Biostrings)
library (readr)
library (dplyr)

#Load the sequences data analyzed in Python
all_sequences <- read_csv(
    "data/processed/validated_sequences.csv",
    show_col_types = FALSE
)

all_sequences$sequence_id <- sprintf(
    "seq_%05d",
    seq_len(nrow(all_sequences))
)

sequences <- all_sequences %>%
    filter(valid_dna == TRUE) 

#Define motifs
motifs <- c(
    start_codon_ATG = "ATG",
    stop_codon_TAA = "TAA",
    stop_codon_TAG = "TAG",
    stop_codon_TGA = "TGA",
    TATA_box = "TATA"
)

#Convert DNA strings into Biostrings objects
dna_set <- DNAStringSet(sequences$sequence)

#Load python results
python_results <- read_csv(
    "results/motif_summary.csv",
    show_col_types = FALSE
)

print(names(python_results))
print(head(python_results))

#Count motifs using Biostrings
r_results <-lapply(names(motifs), function(name) {
    counts <- vcountPattern(
        motifs[[name]],
        dna_set,
        fixed=TRUE
    )
    data.frame(
        sequence_id = sequences$sequence_id,
        motif_name = name,
        r_count = counts
    )
})

r_results <- bind_rows(r_results) 

#Compare results
comparison <- python_results %>%
  dplyr::select(sequence_id, motif_name, count) %>%
  dplyr::left_join(
    r_results,
    by = c("sequence_id", "motif_name")
  ) %>%
  dplyr::mutate(matches = count == r_count)

# Check results
print(table(comparison$matches, useNA = "ifany"))

if (any(is.na(comparison$matches)) ||
    !all(comparison$matches)) {
  stop("Python and R motif counts do not fully agree.")
}

cat("\nSUCCESS: All Python and R motif counts match!\n")