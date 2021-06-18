#Extracts a given region from an alignment
#First genome position is position 1 so use 1 based numbers
#To run python3 extract_alignment_region.py -a alignment.fasta -p1 start_position -p2 end_position -o output.fasta

import argparse
from Bio import AlignIO

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", help = "Input fasta alignment")
    parser.add_argument("-p1", help = "Start of the region to be extracted. 1 based so " + 
                                    "position 1 in the alignment should be given as 1")
    parser.add_argument("-p2", help = "End of the region to be extracted. 1 based so " + 
                                    "position 100 in the alignment should be given as 100")
    parser.add_argument("-o", help = "Output fasta alignment containing extracted region")
    args = parser.parse_args()

    #Import the alignment
    alignment = AlignIO.read(args.a, "fasta")

    outFile = open(args.o, "w")

    #Positions to be extracted
    p1 = int(args.p1) - 1
    p2 = int(args.p2)

    #Iterate through the sequences and write the region to be extracted
    for s in alignment:
        outFile.write(">" + s.id + "\n" + str(s.seq)[p1:p2] + "\n")

    outFile.close()