# sPRuce
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
```

```
bash /Data/corey/Software/sPRuce/sPRuce.sh -q <query_fasta> -t <n_threads> [-s <tree_size>] [-o <outgroup>] [-m <tree_model>] -p prefix
```
[optional parameters]
## Options
- q <query_fasta>: Path to the query fasta file.
- t <n_threads>: Number of threads for parallel processing.
- s <tree_size>: Tree size. Options: basic, top50, decent, large, all, or focus. [Default: basic]
- o <outgroup>: Outgroup. Choose a taxon or allow sPRUce to choose one. [Default: auto].
- m <tree_model>: Tree model for IQ-TREE. [Default: MFP]
- p <prefix>: Prefix for output files (Default: output)"

### tree_size
tree_size might take a little trial and error. These size classifications are relative to size of the reference taxa so even "basic" can produce a tree with hundreds of lineages. 
 - "basic" will include one representive sequence for each family in the same subdivision as the query sequence
 - "top50" uses the top 50 blast hits however there is a qcovhsp cutoff so there may be anywhere from 0-50
 - "decent" will include one representive sequence from each genus in the same subdivision as the query sequence
 - "large" will include one representive sequence from each genus in the same subdivision as the query sequence
 - "all" will include all sequences in the same subdivision as the query sequence.
 - "focus" will include all sequences in the same order as the query sequence.

## Examples
```
sPRuce -q unknown_SSU.fasta -t 6 -prefix CoHo16
sPRuce -q unknown_dinoflagellate.fasta -t 6 -s top50 -o Stramenopiles -prefix CoHo17
```

## Outputs
sPRuce produces compiled PR2 data, alignment files, IQTree files, and a rendered phylogenetic tree pdf with coloured labels. 
Before IQtree, it will also produce a .fasttree so you don't have to wait.

Example of final tree using tree_size = focus
[deor_sPRuce_tree.pdf](https://github.com/coreyholt/sPRuce/files/11631870/deor_sPRuce_tree.pdf)

## Caveats
- Although sPRuce can handle multiple input sequences, it may struggle to choose an appropriate outgroup if they are not _closely related_ (within the same supergroup).
- Some functions rely on random choices so you will likely get a different tree if you rerun sPRuce on the same data. This can be useful if it chooses a weird outgroup. 

## Let me know if you run into any issues!
