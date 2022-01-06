#Plots a presence absence heatmap for a given set of genes given a tree
#To run: RScript genePresenceAbsencePhylogeny.R -t tree -g gene_presence_absence.Rtab -s GenesToPlot -o OutputFile

library(argparse)
library(ggtree)
library(ggplot2)

parser <- ArgumentParser()
parser$add_argument("-t", help = "Newick tree against which genes will be plotted")
parser$add_argument("-g", help = "gene_presence_absence.Rtab from Panaroo")
parser$add_argument("-s", help = "File containing names of genes to be plotted, one per line with no header. These should match the gene names in gene_presence_absence.Rtab")
parser$add_argument("-f", help = "Font size for gene names in heatmap, default 2. This can be adjusted if the names are too large", default = "2")
parser$add_argument("-o", help = "Name of output PDF containing heatmap")
args <- parser$parse_args()

#Import the tree
tree <- read.tree(args$t)
#Import gene_presence_absence.Rtab
genePresenceAbsence <- read.table(args$g,sep="\t",header=TRUE,comment.char="!",check.names=FALSE)
#List of the gene names to be plotted
genesToPlot <- read.table(args$s, sep = "\t", header = FALSE)
outFile <- args$o

#Extract the genes to be plotted
genePresenceAbsencePlot <- genePresenceAbsence[which(genePresenceAbsence[,1] %in% genesToPlot[,1]),]

#Will be filled with the gene presence absence in each sample
geneSummary <- data.frame(matrix(0,ncol=length(genePresenceAbsencePlot[,1]),nrow=(length(genePresenceAbsencePlot[1,])-1)))
rownames(geneSummary) <- names(genePresenceAbsencePlot)[2:length(names(genePresenceAbsencePlot))]
names(geneSummary) <- c(as.matrix(genePresenceAbsencePlot$Gene))

#Iterate through the genes to be plotted and add their presence absence to geneSummary
for (gene in 1:(length(genePresenceAbsencePlot[1,])-1)) {
  geneSummary[gene,] <- genePresenceAbsencePlot[,(gene+1)]}

p <- ggtree(tree)

pdf(outFile)
phylogenyPlot <- gheatmap(p, geneSummary, font.size = as.numeric(args$f), colnames_angle = 45) + scale_fill_continuous(low = "grey", high = "dodgerblue")
print(phylogenyPlot)
dev.off()
