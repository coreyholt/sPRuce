# Install required packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("pr2database", quietly = TRUE)) {
devtools::install_github("pr2database/pr2database")
}
# Check if tidyverse package is installed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

# Check if Biostrings package is installed
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  install.packages("Biostrings")
}

# Load necessary libraries
suppressPackageStartupMessages(library(pr2database))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Biostrings))

pr2 <- pr2_database()

# Filter pr2 database
pr2_ref <- pr2 %>% 
  filter(organelle == "nucleus")%>% 
  #filter(reference_sequence == 1) %>% 
  select(subdivision, 
         order, 
         species, 
         genbank_accession, 
         sequence) %>% 
  distinct(genbank_accession, .keep_all=TRUE) %>% 
  filter(genbank_accession != "FJ228701") %>% 
  mutate_all(~ str_replace_all(., "_", "-"))

# Write as fasta
pr2_ref_seq <- Biostrings::DNAStringSet(pr2_ref$sequence)

# Change sequence headers
names(pr2_ref_seq) <- paste(pr2_ref$subdivision,
                            pr2_ref$order,
                            pr2_ref$genbank_accession,
                            sep="_")

# Save fasta
Biostrings::writeXStringSet(pr2_ref_seq, 
                            "sPRuce_fullDB.fasta",
                            width = 80)