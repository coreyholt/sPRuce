pdf(NULL)

# Check if castor package is installed
if (!requireNamespace("castor", quietly = TRUE)) {
  install.packages("castor")
}
# Check if tidyverse package is installed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}
# Check if RColorBrewer package is installed
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
# Check if BiocManager package is installed
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
BiocManager::install("treeio")
BiocManager::install("ggtree")
}

# Load requried libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(castor))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(RColorBrewer))

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) != 1) {
  stop("Incorrect number of arguments. Usage: Rscript your_script.R input.fasta tree.size outgroup")
}

# Assign input filename and outgroup value
prefix <- args[1]

# Define "type" palette
pal.type <- c(outgroup = "black",
              reference = "black",
              query = "deepskyblue3")

# Load sPRuce tree
tree <- read.iqtree(paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_alignment_trimal.fasta.contree", sep=""))
tree <- as.phylo(tree)

# Load sPRuce metadata
sPRuce.data <- read_csv(paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_metadata.csv", sep=""), show_col_types = FALSE)

# Replace characters that ggtree renames
sPRuce.data$seq_name <- gsub("[():]", "_", sPRuce.data$seq_name)

# Reorder metadata
sPRuce.data.f <- sPRuce.data %>% 
  select(seq_name, subdivision, order, type)

# Isolate outgroup taxa
outgroup <- sPRuce.data %>% 
  filter(type=="outgroup") %>% 
  pull(seq_name)
  
# Find most recent common ancestor of outgroup taxa
outgroup.node <- get_mrca_of_set(tree, outgroup)

# Reroot according to mrca node
tree.root<-root(tree, node=outgroup.node)

# Make ggtree
p1 <- ggtree(tree.root, size=0.4)

# Add metadata to colour query tip labels. Add circular node labels to represent full bootstrap support
p2 <- p1 %<+% sPRuce.data.f  +
  geom_tiplab(mapping=aes(colour=type),
              size=3, offset=0.015)+
  geom_nodelab(aes(subset = !is.na(as.numeric(label)) & as.numeric(label) != 100),
               nudge_x = -0.0125,
               nudge_y = 0.5,
               size=2,
               colour="black")+
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) == 100),
    colour = "black",
    size = 1
  ) +xlim(0,1.1)+
  scale_color_manual(values=pal.type)+ 
  guides(color = "none")

# Define colour palette
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

# If nummber of unique orders is below 16, label orders. Otherwise, label subdivisions
if (length(unique(sPRuce.data.f$order)) < 16) {
  unique_orders <- unique(sPRuce.data.f$order)
  random_colors <- sample(getPalette(length(unique_orders)), replace = FALSE)
  p3 <- p2 +
    geom_tippoint(mapping = aes(x = x + 0.0075, fill = order), 
                  colour = "black", shape = 21, size = 2) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1, 'mm'),
          legend.spacing.y = unit(0.1, 'mm'),
          legend.margin=margin(-0.5,0,0,2, "cm"))+
    scale_fill_manual(values = random_colors)
} else {
  unique_subdivisions <- unique(sPRuce.data.f$subdivision)
  random_colors <- sample(getPalette(length(unique_subdivisions)), replace = FALSE)
  p3 <- p2 +
    geom_tippoint(mapping = aes(x = x + 0.0075, fill = subdivision), 
                  colour = "black", shape = 21, size = 2) +
    theme(plot.margin = margin(1, 1, 1, 1, "cm")) +
    guides(fill = guide_legend(byrow = TRUE)) +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(1, 'mm'),
          legend.spacing.y = unit(0.1, 'mm'),
          legend.margin=margin(-0.5,0,0,2, "cm"))+
    scale_fill_manual(values = random_colors)
}

p4<- p3 + geom_treescale(offset=-1,
                         y=0)

p4
# Save final tree as pdf
ggsave(plot=p4, 
       paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_tree.pdf", sep=""),
       width=25,
       device="pdf",
       height=nrow(sPRuce.data.f)/2.5,
       units = "cm")
