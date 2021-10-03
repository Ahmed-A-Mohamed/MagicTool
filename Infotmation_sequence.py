from Bio import SeqIO
from Bio.SeqUtils import GC

#import os
#mySequence = os.environ.get("seq")
mySequence = input("Enter your sequence in .fasta format!")

# count the reads
myFile = open(f"result.{mySequence}.ods", "w")
count = 0
for rec in SeqIO.parse(mySequence, "fasta"):
    count += 1

myFile.write("Your file contains %i reads" % count)
myFile.write("\n\n")
myFile.close()
##################################################################################
# Details of sequences
# Count letter
# while loop to counter letter
myFile = open(f"result.{mySequence}.ods", "a")
myFile.write(f"")
myFile.write("\n")
kme = input("Do you need Kmer analysis? Y or N")
if kme == "Y" or kme == "y":
    k = input("Enter the number of K-mers ?")

specific = input("Do you need count specific K-mer? Y or N")
if specific == "Y" or specific == "y":
    Pattern = input("Enter the specific K-mer to count")

def FrequentWords(Text, k):
    # your code here
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        if freq[key] == m:
            words.append(key)
    return words


def FrequencyMap(Text, k):
    # your code here.
    freq = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i + k]
        freq[Pattern] = 0
        # hint: your code goes here!
        for j in range(n - k + 1):
            if Text[j:j + k] == Pattern:
                freq[Pattern] = freq[Pattern] + 1
    return freq




for seq_record in SeqIO.parse(mySequence, "fasta"):
    myFile.write(seq_record.description)
    myFile.write("\n")
    myFile.write(repr(seq_record.seq))
    myFile.write("\n")
    myFile.write("The length of this record is ")
    myFile.write(str(len(seq_record))) # Length of the sequence
    myFile.write("\n")
    myFile.write(f"The count of ( A ) => ")
    myFile.write(str((seq_record.seq).count("A")))
    myFile.write("\n")
    myFile.write(f"The count of ( C ) => ")
    myFile.write(str((seq_record.seq).count("C")))
    myFile.write("\n")
    myFile.write(f"The count of ( G ) => ")
    myFile.write(str((seq_record.seq).count("G")))
    myFile.write("\n")
    myFile.write(f"The count of ( T ) => ")
    myFile.write(str((seq_record.seq).count("T")))
    myFile.write("\n")
    myFile.write(f"The count of ( U ) => ")
    myFile.write(str((seq_record.seq).count("U")))
    myFile.write("\n")
    myFile.write("CG content percentage =>  ")
    myFile.write(str(round(GC(seq_record.seq), 2)))
    myFile.write("%")

    if kme == "Y" or kme == "y":
        #k = input("Enter the number of K-mers ?")

        myFile.write("\n")
        myFile.write("The highest K-mer is => ")
        myFile.write(str(FrequentWords(str(seq_record.seq), int(k))))
        myFile.write("\n")
        myFile.write("The all K-mers are => ")
        myFile.write(str(FrequencyMap(str(seq_record.seq), int(k))))

    if specific == "Y" or specific == "y":

        def PatternCount(Text, Pattern):
            count = 0
            for i in range(len(Text) - len(Pattern) + 1):
                if Text[i:i + len(Pattern)] == Pattern:
                    count = count + 1
            return count


       # Pattern = input("Enter the specific K-mer to count")

        myFile.write("\n")
        myFile.write(f"The count of ({Pattern})  => ")
        myFile.write(str(PatternCount(str(seq_record.seq), Pattern)))

    myFile.write("\n\n\n")

myFile.write("\n\n")
myFile.close()
################################################################################
