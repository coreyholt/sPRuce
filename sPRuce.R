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

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 5) {
  stop("Incorrect number of arguments. Usage: Rscript your_script.R input.fasta tree.size outgroup")
}

# Assign input filename and outgroup value
input_blast <- args[1]
tree_size <- args[2]
outgroup <- args[3]
query_fasta <- args[4]
prefix <- args[5]

#prefix <- "CoHo17_sub1"
#input_blast <- paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_queries_DB.nt.blastn", sep="")
#tree_size <- "basic"
#outgroup <- "auto"
#query_fasta <- "CoHo17_Euk18SRRNA_sub1.fasta"

pr2 <- pr2_database() %>% 
  mutate_all(~ str_replace_all(., "_", "-"))

# Read query BLAST results
query.blast <-read_tsv(input_blast, show_col_types = FALSE)

# Retrieve the top matching sequences for each query based on pident and qcovhsp
query_match_vec <- query.blast %>% 
  filter(qcovhsp >= 80) %>% 
  group_by(qseqid) %>% 
  slice_max(pident, n=10) %>% 
  ungroup() %>% 
  mutate(sseqid = str_extract(sseqid, "^[^_]+")) %>%
  distinct(sseqid) %>% 
  pull(sseqid)

#query_match_top <- query_match_vec[1]

# Retrieve the top 20 matching sequences for each query based on pident and qcovhsp
query_match_vec_20 <- query.blast %>% 
  filter(qcovhsp >= 80) %>% 
  group_by(qseqid) %>% 
  slice_max(pident, n=20) %>% 
  ungroup() %>% 
  mutate(sseqid = str_extract(sseqid, "[^_]+$")) %>% 
  distinct(sseqid) %>% 
  pull(sseqid)

# ensure top 20 pr2 hits are in the tree by adding them separately
pr2_ss_top20 <- pr2 %>% 
  filter(genbank_accession %in% query_match_vec_20) %>% 
  filter(organelle == "nucleus" )%>% 
  select(subdivision, 
         order, 
         species, 
         genbank_accession, 
         sequence)

# Retrieve representative PR2 hits that match the query sequences based on subdivision
pr2_ss <- pr2 %>% 
  filter(subdivision %in% query_match_vec) %>% 
  filter(organelle == "nucleus")%>% 
  filter(reference_sequence == 1) %>% 
  select(subdivision, 
         class,
         order,
         family,
         genus,
         species, 
         genbank_accession, 
         sequence)

if (any(query_match_vec == unique(pr2_ss$subdivision))==FALSE){
  pr2_ss <- pr2 %>% 
    filter(subdivision %in% query_match_vec) %>% 
    filter(organelle == "nucleus")%>% 
    select(subdivision, 
           class,
           order,
           family,
           genus,
           species, 
           genbank_accession, 
           sequence)
}

# Generate filtered PR2 hits based on tree size selection
 if (tree_size == "basic"){
  pr2_ss_grouped_filt <- pr2_ss %>% 
    group_by(order, family) %>% 
    group_modify(~sample_n(.x, size = 1)) %>% 
    ungroup() %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence) %>% 
    bind_rows(pr2_ss_top20) %>% 
    distinct(genbank_accession, .keep_all=TRUE)
  
  
} else if (tree_size == "decent") {
  pr2_ss_grouped_filt <- pr2_ss %>% 
    group_by(order, genus) %>% 
    group_modify(~sample_n(.x, size = 1)) %>% 
    ungroup() %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence)%>% 
    bind_rows(pr2_ss_top20) %>% 
    distinct(genbank_accession, .keep_all=TRUE)
  
} else if (tree_size == "large") {
  pr2_ss_grouped_filt <- pr2_ss %>% 
    group_by(order, species) %>% 
    group_modify(~sample_n(.x, size = 1)) %>% 
    ungroup() %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence)%>% 
    bind_rows(pr2_ss_top20) %>% 
    distinct(genbank_accession, .keep_all=TRUE)
  
} else if (tree_size =="all") {
  pr2_ss_grouped_filt <- pr2_ss %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence)%>% 
    bind_rows(pr2_ss_top20) %>% 
    distinct(genbank_accession, .keep_all=TRUE)
  
} else if (tree_size =="top50") {
  top50_accessions <- query.blast %>% 
    group_by(qseqid) %>% 
    filter(qcovhsp >= 80) %>% 
    slice_max(pident, n=50) %>% 
    ungroup() %>% 
    mutate(sseqid = str_extract(sseqid, "[^_]+$")) %>% 
    distinct(sseqid) %>% 
    pull(sseqid)
  
  pr2_ss_grouped_filt <- pr2_ss %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence)%>% 
    filter(genbank_accession %in% top50_accessions) %>% 
    bind_rows(pr2_ss_top20) %>% 
    distinct(genbank_accession, .keep_all=TRUE)
  
} else if (tree_size =="focus") {
  top_accessions <- query.blast %>% 
    group_by(qseqid) %>% 
    filter(qcovhsp >= 80) %>% 
    slice_max(pident, n=1) %>% 
    ungroup() %>% 
    mutate(sseqid = str_extract(sseqid, "[^_]+$")) %>% 
    distinct(sseqid) %>% 
    pull(sseqid)
  
  top_accessions_focus <- pr2 %>% 
    filter(genbank_accession %in% top_accessions) %>% 
    pull(order)
  
  pr2_ss_grouped_filt <- pr2_ss %>% 
    filter(order %in% top_accessions_focus) %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence) %>% 
    bind_rows(pr2_ss_top20) %>% 
    distinct(genbank_accession, .keep_all=TRUE)
}

# Generate outgroup based on outgroup stragey
if (outgroup == "auto") {
  pr2_match_top_division<- pr2 %>% 
    filter(subdivision %in% query_match_vec) %>% 
    distinct(division) %>% 
    pull(division)
  
  if (all(pr2_match_top_division == pr2_match_top_division[1])==FALSE | any(query_match_vec %in% paste(pr2_match_top_division, "-X", sep=""))){
    pr2_match_top_supergroup<- pr2 %>% 
      filter(subdivision %in% query_match_vec) %>% 
      filter(!str_detect(supergroup, "_X")) %>% 
      distinct(supergroup) %>% 
      pull(supergroup)
    
    supergroup_alldivisions <- pr2 %>% 
      filter(supergroup == pr2_match_top_supergroup) %>% 
      distinct(division) %>% 
      filter(!division %in% pr2_match_top_division) %>% 
      pull(division) 
    
    error_file <- paste(prefix,"_sPRuce_output/",prefix,"_error_file.txt", sep="") 
    if (isEmpty(supergroup_alldivisions) == TRUE) {
      cat("Error Code 2", file = error_file)
      stop()
    }
    
    supergroup_alldivisions_sample <- sample(supergroup_alldivisions,1)
    
    pr2_ss.out <- pr2 %>% 
      filter(division == supergroup_alldivisions_sample) %>% 
      filter(reference_sequence == 1) %>% 
      select(subdivision, 
             order, 
             species, 
             genbank_accession, 
             sequence)%>% 
      distinct(species, .keep_all=TRUE) %>% 
      filter(subdivision == sample(unique(subdivision), 1)) %>% 
      sample_n(size = min(8, n()), replace=FALSE) %>% 
      distinct(genbank_accession, .keep_all=TRUE) 
    
    
    } else if (all(pr2_match_top_division == pr2_match_top_division[1])) {
    pr2_ss.out <- pr2 %>% 
    filter(division %in% pr2_match_top_division) %>% 
    filter(!subdivision %in% query_match_vec) %>% 
    filter(reference_sequence == 1) %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence)%>% 
    distinct(species, .keep_all=TRUE) %>% 
    filter(subdivision == sample(unique(subdivision), 1)) %>% 
      sample_n(size = min(8, n()), replace=FALSE) %>% 
      distinct(genbank_accession, .keep_all=TRUE) 
    
}
  }else {
  pr2_ss.out <- pr2 %>% 
    filter(if_any(3:10, ~ str_detect(., outgroup))) %>%
    filter(organelle == "nucleus")%>% 
    filter(reference_sequence == 1) %>% 
    select(subdivision, 
           order, 
           species, 
           genbank_accession, 
           sequence)%>% 
    distinct(species, .keep_all=TRUE) %>% 
    filter(subdivision == sample(unique(subdivision), 1)) %>% 
    sample_n(size = min(8, n()), replace=FALSE) %>% 
    distinct(genbank_accession, .keep_all=TRUE) 
}

# Combine into single data frame
sPRuce.output.df <- rbind(pr2_ss_grouped_filt,
                          pr2_ss.out)

# Write to fasta
sPRuce.output.seq <- Biostrings::DNAStringSet(sPRuce.output.df$sequence)

# Add sequence headers
names(sPRuce.output.seq) <- paste(sPRuce.output.df$subdivision,
                                  sPRuce.output.df$order,
                                  sPRuce.output.df$species,
                                  sPRuce.output.df$genbank_accession,
                           sep="_")

# Save fasta file
Biostrings::writeXStringSet(sPRuce.output.seq, 
                            paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_refs.fasta",sep=""),
                            width = 80)

# Add metadata to generated dataset
sPRuce.output.df <- sPRuce.output.df %>% 
  mutate(seq_name =names(sPRuce.output.seq)) %>% 
  mutate(type = case_when(
    subdivision == unique(pr2_ss.out$subdivision) ~ "outgroup",
    TRUE ~ "reference"
  ))

query_match_vec_1 <- query.blast %>% 
  filter(qcovhsp >= 80) %>% 
  group_by(qseqid) %>% 
  slice_max(pident, n=1) %>% 
  ungroup() %>% 
  mutate(sseqid = str_extract(sseqid, "[^_]+$")) %>% 
  select(qseqid, sseqid) %>% 
  left_join(pr2, by = c("sseqid" = "genbank_accession")) %>% 
    mutate(
      species = NA,
      sequence = NA,
      genbank_accession =NA,
      type="query"
    ) %>% 
  dplyr::rename(seq_name = qseqid) %>% 
  select(subdivision, 
         order, 
         species, 
         genbank_accession, 
         sequence,
         seq_name,
         type) %>% 
  distinct(seq_name, .keep_all = TRUE) 


# Load query fasta file for sequence
fasta <- readDNAStringSet(query_fasta, format="fasta")
names(fasta) <- sapply(strsplit(names(fasta), " "), function(x) x[[1]])
fasta.f <- fasta[query_match_vec_1$seq_name]

# Add query sequences to meta data
query_match_vec_1$sequence <- unlist(sapply(query_match_vec_1$seq_name, function(name) {
  query_seq <- fasta.f[which(names(fasta.f) == name)]
  as.character(query_seq)
}))

# Combine with previous metadata
sPRuce.output.df.query <- rbind(sPRuce.output.df,
                                query_match_vec_1)

# Write as csv file
write_csv(sPRuce.output.df.query,
          paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_metadata.csv", sep=""))
