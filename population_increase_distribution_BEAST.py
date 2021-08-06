#Calculates the date of the first increase in relative genetic diversity in each sampled MCMC step in a BEAST posterior distribution
#Uses the PopSizes and GroupSizes in a BEAST log file to identify if the relative genetic diversity increases
#and the date of the increase if it does
#Takes the .trees and .log files from BEAST
#By default, expects output from BEAST2 in which the log file columns contain PopSize and GroupSize. If using BEAST1, set option -b to 1
#With BEAST1, columns in the log file are expected to contain popSize and groupSize
#To run: python3 population_increase_distribution_BEAST.py -t tree_distribution.trees -l log_file.log -p minimum_percentage_increase -d latest_sample_date -o output_file_name.txt

from io import StringIO
from Bio import Phylo as p
from operator import itemgetter
import re
import argparse

#Removes the header region from a log file
def removeHeader(logFile):
    lines = []

    for line in logFile:
        if (line[0] != "#") and (line != "End;\n") and (line != "end;\n"):
            lines.append(line)
    
    return(lines)

#Extracts the trees from a .trees file
def extractTrees(treesFile):
    lines = []

    for line in treesFile:
        if line[0:4] == "tree":
            lines.append(re.sub(".* ", "", line.strip()))
    
    return(lines)

#Identifies the columns in a log file that correspond to the GroupSizes
def getGroupSizes(logFile, bVersion):
    #Will be filled with the GroupSizes positions in the header
    positions = []

    #Iterate through the column headers and add to positions
    if bVersion == "2":
        for i, column in enumerate(logFile[0].strip().split("\t")):
            if "GroupSizes" in column:
                positions.append(i)
    
    elif bVersion == "1":
        for i, column in enumerate(logFile[0].strip().split("\t")):
            if "groupSize" in column:
                positions.append(i)
    
    return(positions)

#Identifies the columns in a log file that correspond to the PopSizes
def getPopulationSizes(logFile, bVersion):
    #Will be filled with the GroupSizes positions in the header
    positions = []

    #Iterate through the column headers and add to positions
    if bVersion == "2":
        for i, column in enumerate(logFile[0].strip().split("\t")):
            if "PopSizes" in column:
                positions.append(i)
    
    elif bVersion == "1":
        for i, column in enumerate(logFile[0].strip().split("\t")):
            if "popSize" in column:
                positions.append(i)
    
    return(positions)

#Extracts and sorts node heights in a given tree
def getNodeHeights(tree, date):
    #Will be filled with all of the node heights in the tree
    nodeHeights = []

    #Calculate the height of the root node
    rootNodeHeight = float(max(tree.depths().values()))
    
    #Iterate through the clades, calculate their date and append to nodeHeights
    for clade in tree.get_nonterminals():
        nodeHeights.append(date - (rootNodeHeight - float(tree.depths()[clade])))
    
    #Sort the nodes heights into date order
    sortedNodeHeights = sorted(nodeHeights)

    return(sortedNodeHeights)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", help = "The .trees file from BEAST containing the distribution of trees")
    parser.add_argument("-l", help = "The .log file from BEAST")
    parser.add_argument("-p", help = "Minimum percentage increase in relative genetic diversity above baseline " + 
                                    "to define a population increase, default = 100", default = "100")
    parser.add_argument("-d", help = "Date of latest sample as decimal, e.g. 2015.54")
    parser.add_argument("-b", help = "BEAST version used. Can either be 1 or 2, default is 2", default = "2")
    parser.add_argument("-n", help = "Print update every nth tree. Default is 1000", default="1000")
    parser.add_argument("-o", help = "Output file name")
    args = parser.parse_args()

    #Import the log file and remove its header
    log = open(args.l).readlines()
    logFile = removeHeader(log)

    outFile = open(args.o, "w")
    outFile.write("MCMC_state\tIncrease_date\n")

    #Identify the columns in the log file that correspond to the PopSizes and GroupSizes
    groupPositions = getGroupSizes(logFile, args.b)
    populationPositions = getPopulationSizes(logFile, args.b)

    #Remove the header from the log file
    del(logFile[0])

    #Will be incremented with each tree
    logLine = 0

    #Incremented with each tree
    j = 0
    #Incremented with each tree with an increase in relative genetic diversity
    k = 0

    #Iterate through the trees, identify the corresponding log line and determine if and when the relative genetic diversity increased
    with open(args.t) as fileobject:
        for line in fileobject:

            if line[0:4] == "tree":

                #Print update every nth tree
                if j % int(args.n) == 0:
                    print("Analysing tree", j)
                j += 1

                #Read in the line as a tree
                tree = p.read(StringIO(re.sub(".* ", "", line)), "newick")
                #Extract the corresponding line from the log file
                logTree = logFile[logLine]
                #Extract the MCMC state
                MCMCState = logTree.strip().split("\t")[0]

                #Extract the node heights in the tree
                nodeHeight = getNodeHeights(tree, float(args.d))

                #Extract the GroupSizes and PopSizes for the current MCMC step
                groupSizes = [logTree.strip().split("\t")[i] for i in groupPositions][::-1]
                populationSizes = [logTree.strip().split("\t")[i] for i in populationPositions][::-1]
                #The size of the population above which will be considered an increase
                basePopulation = float(populationSizes[0])
                basePopulationIncrease = basePopulation + (basePopulation * (float(args.p)/float(100)))

                #Will change away from None if there is an increase within the current MCMC step
                increaseDate = None

                #Iterate through the adjacent population size pairs, check if the population has increased >= args.p above baseline
                #and identify the number of nodes at the switch point if so
                for populationSize in range(1, len(populationSizes)):
                    if float(populationSizes[populationSize]) > basePopulationIncrease:
                        #Check if the increase occurred in the first group, if it does use the number of nodes in the first group to identify the increase date
                        if populationSize == 1:
                            nodeNumber = int(float(groupSizes[0]))
                        #Sum the nodes up to the increase to identify the increase date
                        else:
                            nodeNumber = sum([int(float(k)) for k in groupSizes[:populationSize]])
                        increaseDate = nodeHeight[int(nodeNumber)]
                        break
                
                if increaseDate:
                    outFile.write(MCMCState + "\t" + str(increaseDate) + "\n")
                    k += 1

                logLine += 1
    
    print("Proportion of trees with an inferred increase in relative genetic diversity of " +
        args.p + "% above baseline (0.0 is none, 1.0 is all trees): " + str(float(k)/float(j)))
    
    outFile.close()