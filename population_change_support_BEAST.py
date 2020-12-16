#Calculates the proportion of sampled MCMC steps that support an increase or decrease in relative genetic diversity within
#a given time window
#Takes the .trees and .log files from BEAST
#Uses the PopSizes and GroupSizes to identify changes in relative genetic diversity and their dates
#Supply the time window of interest with -w and give this 2 arguments - the first is the start of the window of interest and
#the second is the end of the window of interest
#By default, looks for a population increase. To look for a decrease, use --decrease
#Provide the required population increase/decrease with -p. An increase/decrease of at least this level within the window of interest is
#looked for. The relative genetic diversity at the start of the window is used as the baseline and increases/decreases measured from this
#Expects trees and log files from BEAST2 by default. If using BEAST1 output, use -v 1
#To run:
#python3 population_change_support_BEAST.py -t tree_distribution.trees -l log_file.log -p minimum_percentage_increase -d latest_sample_date -w window_start window_end -o output_file_name.txt

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
    parser.add_argument("-w", help = "Time window of interest. Takes 2 decimal number: the start of the window " + 
                                    "and the end of the window, e.g. 2010 2015 will look for a change between " +
                                    "2010 and 2015. To look for a change at any date, set these to dates " + 
                                    "outside the dates covered by the tree",
                                    nargs = 2)
    parser.add_argument("--decrease", help = "Use this option to look for a population decrease between the supplied dates. " + 
                                    "If this option is not supplied, an increase is looked for",
                                    action = "store_true", default = False)
    parser.add_argument("-d", help = "Date of latest sample as decimal, e.g. 2015.54")
    parser.add_argument("-b", help = "BEAST version used. Can either be 1 or 2, default is 2", default = "2")
    parser.add_argument("-n", help = "Print update every nth tree. Default is 1000", default="1000")
    parser.add_argument("-o", help = "Output file name")
    args = parser.parse_args()

    #Import the log file and remove its header
    log = open(args.l).readlines()
    logFile = removeHeader(log)

    #Identify the columns in the log file that correspond to the PopSizes and GroupSizes
    groupPositions = getGroupSizes(logFile, args.b)
    populationPositions = getPopulationSizes(logFile, args.b)

    #Remove the header from the log file
    del(logFile[0])

    #Will be incremented with each tree
    logLine = 0

    #Extract the start and end of the window of interest
    windowStart = float(args.w[0])
    windowEnd = float(args.w[1])

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

                #Will change away from None if there is an increase within the window in the current MCMC step
                changeDate = None
                basePopulation = None

                #Calculate the relative genetic diversity at the start of the window of interest
                for i, groupSize in enumerate(groupSizes):
                    #Check if the end of the current group of nodes is in the window, the first group that is will be the base population
                    if float(nodeHeight[sum(int(float(a)) for a in groupSizes[:(i + 1)])]) >= windowStart:
                        startGroup = i
                        basePopulation = float(populationSizes[i])
                        basePopulationIncrease = basePopulation + (basePopulation * (float(args.p)/float(100)))
                        basePopulationDecrease = basePopulation - (basePopulation * (float(args.p)/float(100)))
                        break
                
                #Check if the tree spans the window
                if basePopulation:
                    #The population changes within the window of interest if its first node is within the window
                    #Iterate through the remaining groups, check if they start in the window, if they do check if they change by the required amount
                    for i, eachGroup in enumerate(groupSizes[(startGroup + 1):]):
                        #The first node in the current window is the sum of the nodes in the previous windows
                        #Check if the switch is within the window of interest
                        if float(nodeHeight[sum(int(float(a)) for a in groupSizes[:(startGroup + i + 1)])]) <= windowEnd:
                            #Check if the group has the required population change
                            if args.decrease:
                                if float(populationSizes[startGroup + i + 1]) < basePopulationDecrease:
                                    changeDate = "Yes"
                            else:
                                if float(populationSizes[startGroup + i + 1]) > basePopulationIncrease:
                                    changeDate = "Yes"
                
                #Check if there was an increase/decrease in the window of interest within this MCMC step
                if changeDate == "Yes":
                    k += 1

                logLine += 1
    
    print("The proportion of trees with a population change in the required window is " + str(float(k)/float(j)))