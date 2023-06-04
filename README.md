                       ___
         ____   _____ /_  )         / \
        |  __ \|  __ \ / /_        /   \
     ___| |__) | |__) |__ _/  __  /_/ \_\
    / __|  ___/|  _  / | | |/ __/ _ \  \_\
    \__ \ |    | | \ \ |_| | (_|  __/ \_\
    |___/_|    |_|  \_\__,_|\___\___|_|  

sPRuce takes a fasta file of 18S rRNA gene sequences and creates a phylogenetic tree based on PR2 data. 

To use sPRuce on jezero
```
conda activate sPRuce
bash /Data/corey/Software/sPRuce/sPRuce.sh -q <query_fasta> -t <n_threads> [-s <tree_size>] [-o <outgroup>] [-m <tree_model>] -p prefix
```
[optional parameters]
## Options
```
 - q <query_fasta>: Path to the query fasta file.
 - t <n_threads>: Number of threads for parallel processing.
 - s <tree_size>: Tree size. Options: basic, top50, decent, large, all, or focus. [Default: basic]
 - o <outgroup>: Outgroup. Choose a taxon or allow sPRUce to choose one. [Default: auto].
 - m <tree_model>: Tree model for IQ-TREE. [Default: MFP]
 - p <prefix>: Prefix for output files
```
### tree_size
tree_size might take a little trial and error. These size classifications are relative to size of the reference taxa so even "basic" can produce a tree with hundreds of lineages. **I'd recommend starting with basic.**
 - "basic" will include one representive sequence for each family in the same subdivision as the query sequence
 - "top50" uses the top 50 blast hits... however there is a qcovhsp cutoff so there may be anywhere from 0-50 (see Example 2)
 - "decent" will include one representive sequence from each genus in the same subdivision as the query sequence
 - "large" will include one representive sequence from each species in the same subdivision as the query sequence
 - "all" will include all sequences in the same subdivision as the query sequence.
 - "focus" will include all sequences in the same order as the query sequence.

**All trees will include the top 20 blast hits from PR2.**

## Examples
```
sPRuce -q one_unknown.fasta -t 6 -prefix CoHo16
sPRuce -q two_unknown_SAR.fasta -t 6 -s top50 -o Stramenopiles -prefix CoHo17
```

## Outputs
sPRuce produces compiled PR2 data, alignment files, IQTree files, and a rendered phylogenetic tree pdf with coloured labels. 
Before IQtree, it will also produce a .fasttree so you don't have to wait.

Tip labes are made up of subdivision_order_species_genbankaccession. 
In rendered pdf trees, lineages will be coloured by order providing there are fewer than 16. If not, they will be coloured by subdivision. 
Node support values are indicated by coloured circles to avoid messy positions. Full support = black, 90-99 = dark grey, 70-89 = light grey, less than 70 = white. 

### Example 1 tree using default tree_size (basic) and outgroup (auto) 
![deor_basic_sPRuce_tree](https://github.com/coreyholt/sPRuce/assets/75506746/a0662852-3aa3-4313-8b38-a90d143841e0)

### Example 2 tree using two different input sequences, tree_size = top50, and -o Stramenopiles
![both_top50_sPRuce_tree](https://github.com/coreyholt/sPRuce/assets/75506746/63fd2c02-a715-4c1a-bdad-b6835c5c84ff)

If you are unhappy with the colours you can rerun the tree script 
``` 
RScript /Data/corey/Software/sPRuce/sPRuce_tree.R <prefix>
```
(using the original prefix)

## Caveats
- sPRuce uses a minimum qcovhsp value of 80 when sorting BLAST hits against PR2.
- Trees are based on PR2 "reference sequences" – a subset of the longest sequences from each group. However, all PR2 sequences are used for BLAST. 
- Although sPRuce can handle multiple input sequences, it may struggle to choose an appropriate outgroup if they are not _closely related_ (within the same supergroup).
- Some functions rely on random choices so you will likely get a different tree if you rerun sPRuce on the same data. This can be useful if it chooses a weird outgroup. 
- The final plot width will scale according to branch length distribution but may cut off very long branching taxa. In which case, check tree files. 

## Let me know if you run into any issues!
