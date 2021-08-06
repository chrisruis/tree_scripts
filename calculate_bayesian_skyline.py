#Extracts relative genetic diversity through time for each sample step in a BEAST posterior distribution
#Divides the tree into specified windows and calculates the relative genetic diversity in each window in each sampled step in the posterior
#The start and end of the examined period are provided with -d1 and -d2, respectively, and -a specifies the number of windows in that range. Therefore
#using -d1 1900 -d2 2000 -a 100 will split time in year long windows
#To run: python calculate_bayesian_skyline.py -l BEAST.log -t BEAST.trees -s LatestSampleDate -d1 EarliestDate -d2 LatestDate -a NumberOfIntervals -o OutputFile

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

#Takes 2 dates and a required number of intervals and returns the start and end of each interval
def getStartEnd(date1,date2,numberIntervals):
    #Will be filled with the start and end of each interval
    intervals = []

    for sampleDate in range(int(numberIntervals)):
        intervals.append([(float(date1)+((float(date2)-float(date1))/float(numberIntervals))*sampleDate),(float(date1)+((float(date2)-float(date1))/float(numberIntervals))*(sampleDate+1))])

    return(intervals)

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
    parser.add_argument("-l", help = "Log file from BEAST")
    parser.add_argument("-t", help = "Trees file from BEAST")
    parser.add_argument("-s", help = "Date of the latest sample")
    parser.add_argument("-d1", help = "The earliest date to be examined")
    parser.add_argument("-d2", help = "The latest date to be examined")
    parser.add_argument("-b", help = "BEAST version used. Can either be 1 or 2, default is 2", default = "2")
    parser.add_argument("-a", help = "The number of windows to be examined, default 100", default = "100")
    parser.add_argument("-n", help = "Print update every nth tree. Default is 1000", default = "1000")
    parser.add_argument("-o", help = "Output file")
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

    outFile = open(args.o,"w")

    #Extract start and end of each interval to be examined
    populationIntervals = getStartEnd(args.d1, args.d2, args.a)
    
    j = 0

    outFile.write("Sample\t" + "\t".join([str(m[0]) for m in populationIntervals]) + "\n")

    #Iterate through the trees, identify the corresponding log line and extract the relative genetic diversity in each window
    with open(args.t) as fileobject:
        for line in fileobject:

            if line[0:4] == "tree":

                #Print update every nth tree
                if j %int(args.n) == 0:
                    print("Tree", j)
                j += 1

                #Read in the line as a tree
                tree = p.read(StringIO(re.sub(".* ", "", line)), "newick")
                #Extract the corresponding line from the log file
                logTree = logFile[logLine]

                #Extract the node heights in the tree
                nodeHeight = getNodeHeights(tree, float(args.s))

                #Extract the GroupSizes and PopSizes for the current MCMC step
                groupSizes = [int(logTree.strip().split("\t")[i]) for i in groupPositions][::-1]
                populationSizes = [logTree.strip().split("\t")[i] for i in populationPositions][::-1]

                #The dates at which the relative genetic diversity changes, starts with the root date
                populationChanges = [float(args.s) - max(tree.depths().values())]

                #Iterate through the nodes where the relative genetic diversity changes and add the dates to populationChanges
                for node in range(len(groupSizes)):
                    populationChanges.append(nodeHeight[(sum(groupSizes[:(node+1)])-1)])
                
                #Extract intervals that correspond to each window of relative genetic diversity in this tree
                populationSamples = []
                for sampleDate in range(len(populationChanges)-1):
                    populationSamples.append([populationChanges[sampleDate], populationChanges[sampleDate+1]])
                
                #Extract the relative genetic diversity for each specified window in this tree
                populationSize = [[] for i in range(int(args.a))]
                for k, sampleInterval in enumerate(populationIntervals):
                    for l,samplePopulation in enumerate(populationSamples):
                         #Check if the relative genetic diversity changes within the current interval
                        if (sampleInterval[0] >= samplePopulation[0]) and (sampleInterval[0] <= samplePopulation[1]) and (sampleInterval[1] >= samplePopulation[1]):
                            if l != (len(populationSamples)-1):
                                #Assign to the relative genetic diversity of the later window
                                populationSize[k] = populationSizes[l+1]
                            else:
                                #If the interval crosses the end of the last window, use the population size value in the last window
                                populationSize[k] = populationSizes[l]
                        #Check if the interval is within the population interval
                        elif sampleInterval[0] >= samplePopulation[0] and sampleInterval[1] <= samplePopulation[1]:
                            populationSize[k] = populationSizes[l]
            
                #Add zeros where the date is not spanned by the tree
                for dateInterval in range(len(populationSize)):
                    if populationSize[dateInterval] == []:
                        populationSize[dateInterval] = "0"
                
                #Write the relative genetic diversity in each window in this tree
                outFile.write("Sample" + str(j) + "\t" + "\t".join(populationSize) + "\n")

                logLine += 1
    
    outFile.close()