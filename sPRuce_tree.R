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
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install treeio package
if (!requireNamespace("treeio", quietly = TRUE)) {
  BiocManager::install("treeio")
}

# Install ggtree package
if (!requireNamespace("ggtree", quietly = TRUE)) {
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
              query = "white")

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
  ) 

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
                            fill=order, 
                            alpha=type),
                size=3,
                geom="label",
                fontface = "bold",
                label.size = NA,
                label.r = unit(0.5, "lines"),
                offset=offset_value)+
    scale_color_manual(values=pal.type)+ 
    scale_alpha_manual(values=c(0.5,1,0.5))+
    guides(color = "none",
           fill = guide_legend(byrow = TRUE,
                               override.aes = list(size = 4,
                                                   color=NA,
                                                   alpha=NA)),
           alpha = "none")+
             theme(plot.margin = margin(1, 1, 1, 1, "cm"),
                   legend.key=element_rect(colour=NA, fill=NA),
                   legend.position = "bottom",
                   legend.spacing.x = unit(1, 'mm'),
                   legend.spacing.y = unit(0.1, 'mm'),
                   legend.margin=margin(-0.5,0,0,2, "cm"))+
    scale_fill_manual(values = random_colors, 
                      na.value = "grey15", 
                      limits = na.omit(unique_orders))
} else {
  unique_subdivisions <- unique(sPRuce.data.f$subdivision)
  random_colors <- sample(getPalette(length(unique_subdivisions)), replace = FALSE)
  p3 <- p2  +
    geom_tiplab(mapping=aes(colour=type,
                            fill=subdivision,
                            alpha=type),
                size=3,
                geom="label",
                fontface = "bold",
                label.size = NA,
                label.r = unit(0.5, "lines"),
                offset=offset_value)+
    scale_color_manual(values=pal.type)+ 
    scale_alpha_manual(values=c(0.5,1,0.5))+
    guides(color = "none",
           fill = guide_legend(byrow = TRUE,
                               override.aes = list(size = 4,
                                                   color=NA,
                                                   alpha=NA)),
           alpha = "none")+
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),
          legend.key=element_rect(colour=NA, fill=NA),
          legend.position = "bottom",
          legend.spacing.x = unit(1, 'mm'),
          legend.spacing.y = unit(0.1, 'mm'),
          legend.margin=margin(-0.5,0,0,2, "cm"))+
    scale_fill_manual(values = random_colors, 
                      na.value = "grey15", 
                      limits = na.omit(unique_subdivisions))
}
p3

p4<- p3 + geom_treescale(y=-1)

p4 

# A (mostly garbage) way to fit everything on while considering long branch lengths
label_coords <- ggplot_build(p4)$data[[7]][c("x", "y")]
max_x <- max(label_coords$x)

p5 <- p4 + xlim(0,max_x*2)
p5

# Save final tree as pdf
ggsave(plot=p5, 
       paste(prefix,"_sPRuce_output/",prefix,"_sPRuce_tree.pdf", sep=""),
       width=25,
       device="pdf",
       height=num_lineages / 2,
       dpi="screen",
       units = "cm",
       limitsize=FALSE)

#prefix<-"telonemaTest"
