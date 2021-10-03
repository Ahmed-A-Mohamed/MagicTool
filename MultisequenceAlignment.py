from Bio.Align.Applications import ClustalwCommandline
import os
mySequence = os.environ.get("seq")

#mySequence = input("Enter your sequence in .fasta format!")

open(mySequence)
cline = ClustalwCommandline("clustalw", infile= mySequence)
cline()