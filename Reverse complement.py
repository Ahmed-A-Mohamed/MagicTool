print("Please confirm that your sequence folder in the same file where this script inside it too.")

from Bio import SeqIO
from Bio.Seq import reverse_complement

#mySequence = input("Enter your sequence in .fasta format!")
import os
mySequence = os.environ.get("seq")

mySeq = open(f"./reverse_complement.{mySequence}.txt", "w")

for seq_record in SeqIO.parse(mySequence, "fasta"):
    mySeq.write(seq_record.name)
    mySeq.write("\n")
    mySeq.write(str(reverse_complement(seq_record.seq)))
    mySeq.write("\n\n")
mySeq.close()