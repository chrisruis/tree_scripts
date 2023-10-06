#Extracts a given region or set of sites from an alignment
#To extract a region, use -p1 and -p2 to specify the start and end of the region to be extracted
#To extract a set of positions, use -p to specify a set of positions or -f to use a file containing positions to be extracted
#The file provided with -f should be 1 position per line with no header
#First genome position is position 1 so use 1 based numbers
#To extract a region from an alignment: python3 extract_alignment_sites.py -a alignment.fasta -p1 start_position -p2 end_position -o output.fasta
#To extract a set of positions: python3 extract_alignment_sites.py -a alignment.fasta -p position1 position2 -o output.fasta
#To extract a set of positions from a file: python3 extract_alignment_sites.py -a alignment.fasta -f positions_file.txt -o output.fasta

import argparse
from Bio import SeqIO, AlignIO

#Check arguments
def check_args(args):
    nA = 0
    if args.p1 and args.p2:
        nA += 1
    if args.p:
        nA += 1
    if args.f:
        nA += 1
    if args.variable:
        nA += 1
    if nA != 1:
        raise RuntimeError("Specify an alignment region to be extracted with -p1 and -p2, or a set of sites to be extracted with -p, " + 
                           "or a file containing sites to be extracted with -f, or --variable to extract variable sites")

#Extract a region from an alignment
def extractAlignmentRegion(align, p1, p2, out):
    #Positions to be extracted
    p1 = int(args.p1) - 1
    p2 = int(args.p2)

    outFile = open(out, "w")

    #Iterate through the sequences and write the region to be extracted
    for r in SeqIO.parse(align, "fasta"):
        outFile.write(">" + r.id + "\n" + str(r.seq)[p1:p2] + "\n")
    
    outFile.close()

#Extract variable sites from an alignment
def extractVariable(alignFile, out, vOut):
    #Import the alignment
    align = AlignIO.read(alignFile, "fasta")

    outV = open(vOut, "w")
    outV.write("Variable_alignment_site,Original_alignment_site,Residues\n")

    #Iterator for variable sites
    v = 1

    #Identify variable sites
    vS = []
    for eS in range(len(align[0].seq)):
        sSet = set(align[:, eS].replace("X", "").replace("-", ""))
        if len(sSet) > 1:
            vS.append(eS)
            outV.write(str(v) + "," + str(eS + 1) + "," + "|".join(sSet) + "\n")
            v += 1
    
    outV.close()
    
    outFile = open(out, "w")

    #Write variable sites
    for eachSeq in align:
        outFile.write(">" + eachSeq.id + "\n")
        for eS in vS:
            outFile.write(eachSeq[eS])
        outFile.write("\n")
    
    outFile.close()

#Extract specified sites from an alignment
def extractSites(align, p, sFile, out):
    #Sites to be extracted
    s = []
    if p:
        for eS in p:
            s.append(int(eS) - 1)
    else:
        with open(sFile) as f:
            for l in f:
                s.append(int(l.strip()) - 1)
        
    outFile = open(out, "w")
    
    #Iterate through the sequences and write the sites to be extracted
    for r in SeqIO.parse(align, "fasta"):
        outFile.write(">" + r.id + "\n")
        for eS in s:
            outFile.write(r.seq[eS])
        outFile.write("\n")
    
    outFile.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help = "Input fasta alignment")
    parser.add_argument("-p1", help = "Start of the region to be extracted. 1 based so " + 
                                    "position 1 in the alignment should be given as 1")
    parser.add_argument("-p2", help = "End of the region to be extracted. 1 based so " + 
                                    "position 100 in the alignment should be given as 100")
    parser.add_argument("-p", help = "Set of positions to be extracted from the alignment. 1 based so " + 
                                    "position 100 in the alignment should be given as 100", nargs = "+")
    parser.add_argument("-f", help = "File containing set of positions to be extracted from the alignment. 1 based so " + 
                                    "position 100 in the alignment should be given as 100")
    parser.add_argument("--variable", help = "Specify to extract variable sites from the alignment. Two output files will be written: " + 
                        "The alignment of variable sites will be written to the file specified with -o. A conversion from the output alignment " + 
                        "position to the input alignment position will also be saved to default file name position_conversion.csv. This " + 
                        "file name can be updated with -vf",
                        action = "store_true",
                        default = False)
    parser.add_argument("-vf", help = "File name to which conversion from variable site alignment to original alignment will be written " + 
                        "if --variable if specified, default position_conversion.csv", default = "position_conversion.csv")
    parser.add_argument("-o", help = "Output fasta alignment containing extracted region")
    args = parser.parse_args()

    #Check arguments
    check_args(args)

    #If a region is to be extracted, extract the region
    if args.p1:
        extractAlignmentRegion(args.a, args.p1, args.p2, args.o)
    #If variable if specified, extract variable sites
    elif args.variable:
        extractVariable(args.a, args.o, args.vf)
    #If sites are specified, extract the sites
    else:
        extractSites(args.a, args.p, args.f, args.o)