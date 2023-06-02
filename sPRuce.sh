#!/bin/bash

# ASCII art
cat << "EOF"
                       ___
         ____   _____ /_  )         / \
        |  __ \|  __ \ / /_        /   \
     ___| |__) | |__) |__ _/  __  /_/ \_\
    / __|  ___/|  _  / | | |/ __/ _ \  \_\
    \__ \ |    | | \ \ |_| | (_|  __/ \_\
    |___/_|    |_|  \_\__,_|\___\___|_|  

EOF

# Default prefix
default_tree_size="basic"
default_outgroup="auto"
default_tree_model="MFP"

# Set default prefix value
tree_size=$default_tree_size
outgroup=$default_outgroup
tree_model=$default_tree_model

# Function to display usage and exit
usage() {
    echo "Usage: $0 -q <query_fasta> -t <n_threads> [-s <tree_size>] [-o <outgroup>] [-m <tree_model>] [-p <prefix>]"
    echo "Options:"
    echo "  -q <query_fasta>       Path to the query fasta file"
    echo "  -t <n_threads>         Number of threads for parallel processing"
    echo "  -s <tree_size>         Tree size: basic, top50, decent, large, or all"
    echo "                         (default: $default_tree_size)"
    echo "  -o <outgroup>          Outgroup: choose taxon of allow sPRuce to choose one"
    echo "                         (default: $default_outgroup)"
    echo "  -m <tree_model>        Tree model for IQ-TREE (default: $default_tree_model)"
    echo "  -p <prefix>            Prefix for output files"
    exit 1
}

# Parse command-line options
while getopts ":q:t:s:o:m:p:" opt; do
    case $opt in
    q)
        query_fasta=$OPTARG
        ;;
    t)
        n_threads=$OPTARG
        ;;
    s)
        tree_size=$OPTARG
        ;;
    o)
        outgroup=$OPTARG
        ;;
    m)
        tree_model=$OPTARG
        ;;
    p)
        prefix=$OPTARG
        ;;
    \?)
        echo "Invalid option: -$OPTARG"
        usage
        ;;
    :)
        echo "Option -$OPTARG requires an argument."
        usage
        ;;
    esac
done

# Check required options
if [[ -z $query_fasta || -z $n_threads || -z $prefix ]]; then
    echo "Error: Missing required options."
    usage
fi

# Create output directory
output_dir="${prefix}_sPRuce_output"
mkdir -p "$output_dir"

# Check if sPRuce_fullDB.fasta doesn't exist, and if so, run sPRUce_DB.R
if [ ! -f "sPRuce_fullDB.fasta" ]; then
    echo "Building sPRuce database..."
    Rscript /Users/coreyholt/Bioinformatics/GitHub/sPRuce/sPRuce_DB.R >/dev/null 2>&1
fi

# Move sPRuce_fullDB.fasta to the output directory
mv sPRuce_fullDB.fasta "$output_dir/"

# Create blast database from sPRuce_fullDB.fasta
makeblastdb -in "$output_dir/sPRuce_fullDB.fasta" -dbtype 'nucl' -parse_seqids -out "$output_dir/sPRuce_DB" >/dev/null 2>&1

# Run blastn with query file and save the output in blastn format
echo "Running blastn..."
blastn -query "$query_fasta" -db "$output_dir/sPRuce_DB" -outfmt '6 qseqid sseqid pident qcovhsp' -num_threads "$n_threads" -evalue 1e-5 > "$output_dir/${prefix}_sPRuce_queries_DB.nt.blastn"

# Add header to the blast output
sed '1s/^/qseqid\tsseqid\tpident\tqcovhsp\n/' "$output_dir/${prefix}_sPRuce_queries_DB.nt.blastn" > "$output_dir/${prefix}_sPRuce_queries_DB.nt.blastn.header"
mv "$output_dir/${prefix}_sPRuce_queries_DB.nt.blastn.header" "$output_dir/${prefix}_sPRuce_queries_DB.nt.blastn"

# Run sPRuce_DB.R with additional arguments
echo "Running sPRuce..."
Rscript /Users/coreyholt/Bioinformatics/GitHub/sPRuce/sPRuce.R "$output_dir/${prefix}_sPRuce_queries_DB.nt.blastn" "$tree_size" "$outgroup" "$query_fasta" "$prefix" >/dev/null 2>&1

# Concatenate sPRuce_refs.fasta and the query file
cat "$output_dir/${prefix}_sPRuce_refs.fasta" "$query_fasta" > "$output_dir/${prefix}_sPRuce_refs_queries.fasta"

# Remove duplicate sequences from the concatenated file
echo "Removing duplicate sequences..."
awk '/^>/{f=!d[$1];d[$1]=1}f' "$output_dir/${prefix}_sPRuce_refs_queries.fasta" > "$output_dir/${prefix}_sPRuce_refs_queries_uniq.fasta"
seqkit rmdup --quiet -s -d "$output_dir/${prefix}_sPRuce_refs_queries_uniq_dups.fasta" < "$output_dir/${prefix}_sPRuce_refs_queries_uniq.fasta" > "$output_dir/${prefix}_sPRuce_refs_queries_uniq2.fasta"

# Check if any query sequences were removed and add them back to the file
while IFS= read -r line
do
  if [[ $line == '>'* ]]; then
    # Extract sequence ID from the header line
    seq_id=$(echo "$line" | sed 's/>//')
    # Search for the sequence in file_b
    grep -q "$seq_id" "$output_dir/${prefix}_sPRuce_refs_queries_uniq_dups.fasta"
    if [ $? -eq 0 ]; then
      # If sequence found in file_b, append it to file_c
      grep -A 1 "$seq_id" "$query_fasta" >> "$output_dir/${prefix}_sPRuce_refs_queries_uniq2.fasta"
      echo "$line" >> "$output_dir/${prefix}_sPRuce_refs_queries_uniq2.fasta"
    fi
  fi
done < "$query_fasta"

# Perform multiple sequence alignment using MAFFT
echo "Performing multiple sequence alignment with MAFFT..."
mafft --quiet --maxiterate 1000 --genafpair --reorder --thread "$n_threads" "$output_dir/${prefix}_sPRuce_refs_queries_uniq2.fasta" > "$output_dir/${prefix}_sPRuce_alignment.fasta"

# Trim the alignment using Trimal
echo "Trimming the alignment with Trimal..."
trimal -in "$output_dir/${prefix}_sPRuce_alignment.fasta" -out "$output_dir/${prefix}_sPRuce_alignment_trimal.fasta" -gt 0.1 -st 0.001 -keepheader >/dev/null 2>&1

# Build a phylogenetic tree using FastTree
echo "Building a phylogenetic tree with FastTree..."
FastTree -quiet -nt -gtr "$output_dir/${prefix}_sPRuce_alignment_trimal.fasta" > "$output_dir/${prefix}_sPRuce_alignment_trimal.fasttree"

# Perform maximum likelihood phylogenetic analysis using IQ-TREE
echo "Performing maximum likelihood phylogenetic analysis with IQ-TREE..."
iqtree -s "$output_dir/${prefix}_sPRuce_alignment_trimal.fasta" -st DNA -m "$tree_model" -bb 1000 -alrt 1000 -nt AUTO -ntmax "$n_threads" -quiet 

rm "$output_dir/sPRuce_DB"* "$output_dir/sPRuce_fullDB.fasta" "$output_dir/${prefix}_sPRuce_refs_queries_uniq2.fasta" "$output_dir/${prefix}_sPRuce_refs_queries.fasta" "$output_dir/${prefix}_sPRuce_refs_queries_uniq_dups.fasta" "$output_dir/${prefix}_sPRuce_refs_queries_uniq.fasta"

# Run sPRuce_tree.R with additional arguments
echo "Making a jazzy pdf..."
Rscript /Users/coreyholt/Bioinformatics/GitHub/sPRuce/sPRuce_tree.R $prefix >/dev/null 2>&1

echo "All Done!"
