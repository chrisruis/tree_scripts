#Calculates an association index for continuous traits
#Takes a tree. The trait can either be in the tree tip labels after the final _ (use --tip_label in this case) or
#can be given as a separate file with 2 or more columns. The first column needs to contain the tip names as they appear in
#the tree. The subsequent columns should contain the traits of interest. The association index will be calculated for each column

from Bio import Phylo
import pandas as pd
import numpy as np
import random
import argparse

#Extracts labels from a tree that are after the last underscore
#Returns a dictionary with tip names as keys and traits as values
def getTreeLabels(tree):
    lDict = dict()
    lDict["Label"] = dict()

    for tip in tree.get_terminals():
        lDict["Label"][tip.name[:tip.name.rindex("_")]] = tip.name.split("_")[-1]

    return(lDict)

#Extracts labels from a given csv file
def getCsvLabels(labels):
    tips = pd.read_csv(labels)

    #Name of the taxon column
    tN = tips.columns[0]

    lDict = dict()
    for c in tips.columns[1:]:
        lDict[c] = dict()

        for tip in range(tips.shape[0]):
            lDict[c][tips[tN][tip]] = tips[c][tip]
    
    return(lDict)

#Calculates variance for a tree and trait
def getTraitVariance(tree, l, trait, tip_label):
    #Total variance
    traitVariance = float(0)

    for clade in tree.get_nonterminals():
        #Do not analyse the root
        if len(tree.get_path(clade)) != 0:
            #Will be filled with the trait in each tip in the current clade
            tipTrait = list()
            for t in clade.get_terminals():
                if tip_label:
                    tipTrait.append(float(l[trait][t.name[:t.name.rindex("_")]]))
                else:
                    tipTrait.append(l[trait][t.name])
            traitVariance += np.var(tipTrait)
    
    return(traitVariance)

#Calculates the continuous association for a tree for a given set of traits
def continuousAI(tree, bootstraps, labels, tip_label):
    #If the labels are in the tree, extract them from the tree
    if tip_label:
        l = getTreeLabels(tree)
    #Import the labels csv file and extract each column
    else:
        l = getCsvLabels(labels)
    
    #Iterate through the traits to test
    #Iterate through the tree and calculate the variance at each internal node, add to traitVariance
    for trait in l:
        cAI = getTraitVariance(tree, l, trait, tip_label)

        #Calculate the bootstrap continuous association index
        bAI = list()
        for b in range(int(bootstraps)):
            #Assign the trait to tips randomly
            tNames = list(l[trait].keys())
            tValues = random.sample(list(l[trait].values()), len(list(l[trait].values())))

            bDict = dict()
            bDict["Label"] = dict()
            for i in range(len(tNames)):
                bDict["Label"][tNames[i]] = tValues[i]
            bAI.append(getTraitVariance(tree, bDict, "Label", tip_label))
        
        #Number of bootstraps with variance at least as small as real data
        nB = 0
        for eB in bAI:
            if eB <= cAI:
                nB += 1
        
        print("Trait:", trait)
        print("Variance association index with real data:", cAI)
        print("Mean variance association index with bootstraps:", sum(bAI)/len(bAI))
        print("Proportion of bootstraps with variance association index at least as small as real data (p-value):", float(nB)/float(bootstraps))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-t",
                        "--tree",
                        dest = "tree",
                        required = True,
                        help = "Newick format phylogenetic tree")
    parser.add_argument("-b",
                        "--bootstraps",
                        dest = "bootstraps",
                        help = "Number of bootstrap samples to be carried out, default = 1000",
                        default = "1000")
    
    label_parser = parser.add_mutually_exclusive_group(required = True)
    label_parser.add_argument("-l",
                            "--labels",
                            dest = "labels",
                            help = "csv file containing traits. The first column needs to contain the tip names " + 
                            "as they appear in the tree. The remaining columns contain traits to be analysed. The " + 
                            "index will be calculated on each column separately. It is necessary to specify either " +
                            "--tip_label or provide a labels file with -l, not both")
    label_parser.add_argument("--tip_label",
                            dest = "tip_label",
                            help = "Specify this if the trait to be analysed is after the last _ in sequence names in " +
                            "the tree. It is necessary to specify either --tip_label or provide a labels file with -l, not both",
                            action = "store_true",
                            default = False)
    
    args = parser.parse_args()

    #Import the tree
    tree = Phylo.read(args.tree, "newick")

    continuousAI(tree, args.bootstraps, args.labels, args.tip_label)