#Extracts a subset of tips from a phylogenetic tree

library(argparse)
library(ape)

parser <- ArgumentParser()
parser$add_argument("-t", help = "Newick tree that will be filtered")
parser$add_argument("-n", help = "Samples to be kept, text file with 1 sample per line and no header")
parser$add_argument("-o", help = "Name of output newick tree")
args <- parser$parse_args()

#Import the tree
tree <- read.tree(args$t)

#Import the tips to be extracted
tipF <- read.csv(args$n, header = FALSE, stringsAsFactors = FALSE)
tips <- tipF[,1]

#Extract the subset of tips from the tree
treeS <- keep.tip(tree, tips)

#Write the tree
write.tree(treeS, file = args$o)