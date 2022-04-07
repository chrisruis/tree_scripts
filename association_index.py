#!/usr/bin/python

#Calculates the association index for a trait of interest on a phylogenetic tree and calculates the statistical significance of that association

description = ""

from Bio import Phylo as p
from collections import Counter
import operator
import random
import statistics
import argparse

def calculateAssociationIndex(phylogeny): #This function takes a phylogeny and calculates the association index of the discrete character at the phylogeny tips
    associationIndex = 0.0 #Will be increased with each internal node that is analysed
    
    for clade in phylogeny.get_nonterminals(): #Iterate through the internal nodes in the phylogeny
        trait = [] #Will be filled with the trait of each tip in the clade
        tipNumber = 0 #Will be increased to the number of tips present in the clade
        for tip in clade.get_terminals(): #Iterate through the tips in the clade
            trait.append(str(tip).split("_")[-1])
            tipNumber += 1
        traitNumber = Counter(trait) #Count the number of occurrences of each trait
        maximumFrequency = max(traitNumber.values()) #Assign to the number of tips with the most frequent trait
        cladeAssociationIndex = (1.0-(float(maximumFrequency)/float(tipNumber)))/((2.0**float(tipNumber))-1) #Calculate the association index of the clade
        associationIndex += cladeAssociationIndex

    return associationIndex

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-t", help = "File path to newick phylogenetic tree with the trait of interest after the last _ in each tip")
    parser.add_argument("-b", help = "Number of bootstraps, default 1000", default="1000")
    args = parser.parse_args()

    phylogeny = p.read(args.t,"newick") #Import the newick phylogeny
    phylogenyAssociationIndex = calculateAssociationIndex(phylogeny)

    phylogenyTipTrait = [] #Will be filled with the trait of each tip in the phylogeny
    for tip in phylogeny.get_terminals(): #Iterate through the tips
        phylogenyTipTrait.append(str(tip).split("_")[-1])

    bootstrapAssociationIndex = [] #Will be filled with the association index for each bootstrap run

    for bootstrap in range(int(args.b)): #Iterate through the bootstrap replicates
        bootstrapSample = random.sample(phylogenyTipTrait,len(phylogenyTipTrait)) #Randomly sample the traits without replacement
        bootstrapPhylogeny = phylogeny
        for i,phylogenyTip in enumerate(bootstrapPhylogeny.get_terminals()): #Iterate through the tips in the bootstrap phylogeny
            phylogenyTip.name = "_" + bootstrapSample[i] #Assign the tip to the ith trait
        bootstrapAssociationIndex.append(calculateAssociationIndex(bootstrapPhylogeny))
    
    numberBootstraps = 0 #Will be increased if the bootstrap has a stronger association index than the real data
    for bootstrapIndex in bootstrapAssociationIndex: #Iterate through the bootstrap association indices
        if bootstrapIndex <= phylogenyAssociationIndex: #Check if the bootstrap association is smaller than the real data
            numberBootstraps += 1
    proportionBootstraps = float(numberBootstraps)/float(args.b) #Calculate the proportion of bootstraps with an association index as strong as the real data

    print("Association Index of the phylogeny = " + str(phylogenyAssociationIndex) + "\nMedian bootstrap Association Index = " + str(statistics.median(bootstrapAssociationIndex)) + "\nP-value on the association = " + str(proportionBootstraps))
