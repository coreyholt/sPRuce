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
              query = "black")

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

# Check if sPRuce.data.f$type is "query" and update subdivision and order accordingly
sPRuce.data.f$subdivision <- ifelse(sPRuce.data.f$type == "query", NA, sPRuce.data.f$subdivision)
sPRuce.data.f$order <- ifelse(sPRuce.data.f$type == "query", NA, sPRuce.data.f$order)

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
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) == 100),
    fill = "black",
    colour="black",
    pch=21,
    size = 2
  ) +
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 100 & as.numeric(label) >= 90),
    fill = "grey25",
    colour="black",
    pch=21,
    size = 2
  ) +
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 90 & as.numeric(label) >= 70),
    fill = "grey70",
    colour="black",
    pch=21,
    size = 2
  ) +
  geom_nodepoint(
    aes(subset = !is.na(as.numeric(label)) & as.numeric(label) < 70),
    fill = "white",
    colour="black",
    pch=21,
    size = 2
  ) +
  coord_cartesian(xlim=c(0,1.1))
p2
# Define colour palette
getPalette = colorRampPalette(brewer.pal(12, "Paired"))

# Calculate the number of lineages in the tree
num_lineages <- length(tree$tip.label)

offset_value <- ifelse(max(tree$edge.length)>7.4, -0.035,
                       ifelse(max(tree$edge.length)>7, -0.025,0))

# If nummber of unique orders is below 16, label orders. Otherwise, label subdivisions
if (length(unique(sPRuce.data.f$order)) < 16) {
  unique_orders <- unique(sPRuce.data.f$order)
  random_colors <- sample(getPalette(length(unique_orders)), replace = FALSE)
  p3 <- p2  +
    geom_tiplab(mapping=aes(colour=type,
                            fill=order),
                size=3,
                geom="label",
                fontface = "bold",
                label.size = 0.5,
                alpha=0.5,
                offset=offset_value)+
    scale_color_manual(values=pal.type)+ 
    guides(color = "none",
           fill = guide_legend(byrow = TRUE))+
             theme(plot.margin = margin(1, 1, 1, 1, "cm"),
                   legend.key=element_rect(colour="black"),
                   legend.position = "bottom",
                   legend.spacing.x = unit(1, 'mm'),
                   legend.spacing.y = unit(0.1, 'mm'),
                   legend.margin=margin(-0.5,0,0,2, "cm"))+
    scale_fill_manual(values = random_colors, 
                      na.value = "white", 
                      limits = na.omit(unique_orders))
} else {
  unique_subdivisions <- unique(sPRuce.data.f$subdivision)
  random_colors <- sample(getPalette(length(unique_subdivisions)), replace = FALSE)
  p3 <- p2  +
    geom_tiplab(mapping=aes(colour=type,
                            fill=subdivision),
                size=3,
                geom="label",
                fontface = "bold",
                label.size = 0.5,
                alpha=0.5,
                offset=offset_value)+
    scale_color_manual(values=pal.type)+ 
    guides(color = "none",
           fill = guide_legend(byrow = TRUE))+
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          legend.key=element_rect(colour="black"),
          legend.position = "bottom",
          legend.spacing.x = unit(1, 'mm'),
          legend.spacing.y = unit(0.1, 'mm'),
          legend.margin=margin(-0.5,0,0,2, "cm"))+
    scale_fill_manual(values = random_colors, 
                      na.value = "white", 
                      limits = na.omit(unique_subdivisions))
}
p3

p4<- p3 + geom_treescale(y=-1)

p4

if (num_lineages > 200){
# Calculate the desired height of the plot based on the number of lineages
plot_height <- num_lineages / 1.75
plot.width <- 24 + (quantile(tree$edge.length,0.75)*50)
} else if (num_lineages < 100){
plot_height <- num_lineages / 1.5
plot.width <- 24 + (quantile(tree$edge.length,0.75)*50)
}
# Save final tree as pdf
ggsave(plot=p4, 
       paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_tree.pdf", sep=""),
       width=plot.width,
       device="pdf",
       height=plot_height,
       dpi="screen",
       units = "cm",
       limitsize=FALSE)

#prefix<-"both_basic"
