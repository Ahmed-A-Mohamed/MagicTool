from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import os
from tkinter import messagebox

window = Tk()
window.title("Magic Tool")
window.geometry('1000x600')

menu = Menu(window)

new_item = Menu(menu)

new_item.add_command(label='About')

new_item.add_separator()

new_item.add_command(label='Help')
new_item.add_separator()

new_item.add_command(label='Exit')

menu.add_cascade(label='File', menu=new_item)
window.config(menu=menu)

tab_control = ttk.Notebook(window)

tab1 = ttk.Frame(tab_control)

tab2 = ttk.Frame(tab_control)
tab3 = ttk.Frame(tab_control)
tab4 = ttk.Frame(tab_control)
tab5 = ttk.Frame(tab_control)
tab6 = ttk.Frame(tab_control)

tab_control.add(tab1, text='Home')

tab_control.add(tab2, text='Pairwise Alignment (DNA)')
tab_control.add(tab3, text='Pairwise Alignment (Protein)')
tab_control.add(tab4, text='Multiple Sequence Alignment')
tab_control.add(tab5, text='Phylogeny')
tab_control.add(tab6, text='Sequence Informaion')


tab_control.pack(expand=1, fill='both')


#############################tab1##########################
frame1 = LabelFrame(tab1, text="Fasta Converter ", padx=200, pady=20)
frame1.pack(padx=1, pady=10)

def search1():
    curr_directory = os.getcwd()
    responce = messagebox.showinfo( "This is information massage", "By clicked okay you wil choose your fastq file and convert it to fasta")
    Label( frame1,text = responce)
    frame1.filename = filedialog.askopenfilename(initialdir= curr_directory, title="add your sequence")
    my_label = Label(frame1, text=frame1.filename)
    my_label.grid()
    return frame1.filename


def convert():
    os.system('''
            cat {seq} | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" > {seq}.fasta 
            '''.format(seq=search1()))

my_btn1 = Button(frame1, text ="Fasta Convertar", command= convert,  borderwidth= 5)
my_btn1.grid()

frame2 = LabelFrame(tab1, text="Download Sequences ", padx=200, pady=20)
frame2.pack(padx=20, pady=20)

txt = Entry(frame2, width=30,  borderwidth= 5)
txt.grid(column=0, row=0)
def clicked():
    res2= messagebox.showinfo("This is information message","By clicked okay you wil download DNA sequence using the accession number in the Entry box")
    os.system('''
    efetch -db nucleotide -id {DNA} -format fasta  > {DNA}.fasta
    '''.format(DNA=txt.get()))
my_btn2= Button(frame2, text= "Download DNA Sequence", command= clicked,  borderwidth= 5)
my_btn2.grid(column=1, row=0)

txt2 = Entry(frame2, width=30, borderwidth= 5)
txt2.grid(column=0, row=2)
def clicked1():
    res3=messagebox.showinfo("This is information message","By clicked okay you wil download Protein sequence using the accession number in the Entry box")
    os.system('''
    efetch -db protein -id {protein} -format fasta  > {protein}.fasta
    '''.format(protein=txt2.get()))
my_btn2= Button(frame2, text= "Download Protein Sequence", command= clicked1,  borderwidth= 5)
my_btn2.grid(column=1, row=2)

txt33 = Entry(frame2, width=30, borderwidth= 5)
txt33.grid(column=0, row=3)
def prot():
    res3=messagebox.showinfo("This is information message","By clicked okay you wil download Protein sequence using the accession number of Nuclotide Sequence in the Entry box")
    os.system('''
    elink -db nuccore  -id {nuclID} -target protein | efetch -format fasta > {nuclID}.protein.fasta
    '''.format(nuclID=txt33.get()))

my_btn2 = Button(frame2, text="Download Protein from Nuclotide ID", command=prot, borderwidth=5)
my_btn2.grid(column=1, row=3)
############################################tab 2 seq to seq#########################################################
frame3 = LabelFrame(tab2, text="Pairwise alignment DNA Sequence to another ", padx=150, pady=20)
frame3.pack(padx=20, pady=20)

def search5():
    curr_directory = os.getcwd()
    res3 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your 2 Files in Fasta formate to make Pairwise alignment ")
    frame3.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame3, text=frame3.filename)
    my_label.grid()
    return frame3.filename

def search51():
    curr_directory = os.getcwd()
    frame3.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame3, text=frame3.filename)
    my_label.grid()
    return frame3.filename

def pairwise():
    os.system('''
     blastn -query {first} -subject {second} -out {first}.PairwiseAlignment.txt 
    '''.format(first=search5(), second=search51()))
mybtn6= Button(frame3, text="Pairwise Alignment 2 DNA Sequences", command =pairwise,  borderwidth= 5)
mybtn6.grid(column=0, row =0)
########################################tab2 sequence to your database

frame4 = LabelFrame(tab2, text="Pairwise alignment DNA Sequence to your Database ", padx=150, pady=20)
frame4.pack(padx=20, pady=20)

def search6():
    curr_directory = os.getcwd()
    res3 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate and make Pairwise alignment against your Database that is in the same location of program and her name in the Entry Box ")
    frame4.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame4, text=frame4.filename)
    my_label.grid()
    return frame4.filename

txt2 = Entry(frame4, width=50,  borderwidth= 5)
txt2.grid(column=0, row=1)
def pairwise2():
    os.system('''
         ./blastn -db {dbname} -query {query} -out {query}.PairwiseAlignment.txt
        '''.format(dbname=txt2.get(), query=search6()))
my_btn6= Button(frame4, text= "Alignment Certain Database ", command =pairwise2,  borderwidth= 5)
my_btn6.grid(column= 1, row = 1)
###########################################tab 2 seq to Human Daabase####################
frame5 = LabelFrame(tab2, text="Pairwise alignment DNA Sequence to Human Database ", padx=150, pady=20)
frame5.pack(padx=20, pady=20)
def search7():
    curr_directory = os.getcwd()
    res4 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate and make Pairwise alignment aginst Human Genoum")
    frame5.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame5, text=frame5.filename)
    my_label.grid()
    return frame5.filename
def pairwise3():
    os.system('''
         ./blastn -db Human_Genoum -query {query2} -out {query2}.PairwiseAlignment.txt
        '''.format(query2=search7()))

my_btn7= Button(frame5, text= "Alignment on Human Genome", command = pairwise3,  borderwidth= 5)
my_btn7.grid(column= 1, row =0 )

def search200():
    curr_directory = os.getcwd()
    res4 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate and make Pairwise alignment aginst Human Genoum")
    frame5.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame5, text=frame5.filename)
    my_label.grid()
    return frame5.filename
txt55 = Entry(frame5, width=30,  borderwidth= 5)
txt55.grid(column=0, row=1)
def pairwisep33():
    os.system('''
         ./blastp -db Human_Protein -query {query22} -outfmt 6 -evalue {e} -out {query22}.pairwise F.txt
        '''.format(query22=search23(), e= txt55.get()))

my_btn8= Button(frame5, text= "Alignment with e-value Filtration", command = pairwisep33,  borderwidth= 5)
my_btn8.grid(column=1 , row =1 )

def indxG():
    res56 = messagebox.askquestion("This is information message", "By clicked okay you wil make Human protein Database to be used as a reference, Please Note that you only need to use this button only one time")
    if res56 == "yes":
        os.system('''
           wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
           gunzip hg38.fa.gz
           makeblastdb -in hg38.fa -dbtype nucl -out Human_Genoum
           ''')
    else:
        Label(frame5, text="process stop").grid()


my_btn6= Button(frame5, text= "Human Genome indexing ", command =indxG,  borderwidth= 5)
my_btn6.grid(column= 1, row = 2)
###################################################tab3 pairwise protein###########################

frame50 = LabelFrame(tab3, text="Pairwise alignment DNA Sequence to another ", padx=150, pady=20)
frame50.pack(padx=20, pady=20)

def search20():
    curr_directory = os.getcwd()
    res3 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your 2 Files in Fasta formate to make Pairwise alignment ")
    frame50.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame50, text=frame50.filename)
    my_label.grid()
    return frame50.filename

def search21():
    curr_directory = os.getcwd()
    frame50.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame50, text=frame50.filename)
    my_label.grid()
    return frame50.filename

def pairwisep():
    os.system('''
     ./blastp -query {first} -subject {second} -out {first}.PairwiseAlignment.txt 
    '''.format(first=search20(), second=search21()))
mybtn20= Button(frame50, text="Pairwise Alignment 2 Protein Sequences", command =pairwisep,  borderwidth= 5)
mybtn20.grid(column=0, row =0)
########################################tab2 sequence to your database

frame51 = LabelFrame(tab3, text="Pairwise alignment Protein Sequence to your Database ", padx=150, pady=20)
frame51.pack(padx=20, pady=20)

def search22():
    curr_directory = os.getcwd()
    res3 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate and make Pairwise alignment against your Database that is in the same location of program and her name in the Entry Box ")
    frame51.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame51, text=frame51.filename)
    my_label.grid()
    return frame51.filename

txt51 = Entry(frame51, width=50,  borderwidth= 5)
txt51.grid(column=0, row=1)
def pairwisep2():
    os.system('''
         ./blastp -db {dbname} -query {query} -out {query}.PairwiseAlignment.txt
        '''.format(dbname=txt51.get(), query=search22()))
my_btn6= Button(frame51, text= "Alignment Protein to Certain Database ", command =pairwisep2, borderwidth= 5)
my_btn6.grid(column= 1, row = 1)
###########################################tab 2 seq to Human Protein Daabase####################
frame52 = LabelFrame(tab3, text="Pairwise alignment Protein Sequence to Human Database ", padx=150, pady=20)
frame52.pack(padx=20, pady=20)
def search23():
    curr_directory = os.getcwd()
    res4 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate and make Pairwise alignment aginst Human Genoum")
    frame52.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame52, text=frame52.filename)
    my_label.grid()
    return frame52.filename
def pairwisep3():
    os.system('''
         ./blastp -db Human_Protein -query {query22} -out {query22}.pairwise.txt
        '''.format(query22=search23()))

my_btn7= Button(frame52, text= "Alignment on Human Protein", command = pairwisep3,  borderwidth= 5)
my_btn7.grid(column= 2, row =0 )
def search27():
    curr_directory = os.getcwd()
    res4 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate and make Pairwise alignment aginst Human Genoum")
    frame52.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame52, text=frame52.filename)
    my_label.grid()
    return frame52.filename
txt53 = Entry(frame52, width=30,  borderwidth= 5)
txt53.grid(column=1, row=1)
def pairwisep33():
    os.system('''
         ./blastp -db Human_Protein -query {query22} -outfmt 6 -evalue {e} -out {query22}.pairwise F.txt
        '''.format(query22=search23(), e= txt53.get()))

my_btn8= Button(frame52, text= "Alignment on Human Protein with Filtration", command = pairwisep33,  borderwidth= 5)
my_btn8.grid(column=2 , row =1 )
def indx():
    res55 = messagebox.askquestion("This is information message", "By clicked okay you wil make Human protein Database to be used as a reference, Please Note that you only need to use this button only one time")
    if res55 == "yes":
        os.system('''
        wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz
        gunzip human.1.protein.faa.gz
        makeblastdb -in human.1.protein.faa -dbtype prot -out Human_Protein
        ''')
    else:
        Label(frame52, text="process stop"). grid()

my_btn6= Button(frame52, text= "Human protein indexing ", command =indx,  borderwidth= 5)
my_btn6.grid(column= 2, row = 2)
#####################################################Tab4 MSA

frame6 = LabelFrame(tab4, text="Multiple Sequence Alignment ", padx=350, pady=20)
frame6.pack(padx=20, pady=20)

def search():
    curr_directory = os.getcwd()
    res4 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate to make Multiple Sequence ALignment")
    frame6.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame6, text=frame6.filename)
    my_label.grid()
    return frame6.filename

def MSA():
    from Bio.Align.Applications import ClustalwCommandline
    mySequence = search()
    open(str(mySequence))
    cline = ClustalwCommandline("clustalw", infile=mySequence)
    cline()
my_btn2= Button(frame6, text="Clustal w", command= MSA, padx=40, pady=20,  borderwidth= 5).grid(row =3, column=2)

def search2():
    curr_directory = os.getcwd()
    res4 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your DNA sequence File in Fasta Formate to make Multiple Sequence ALignment")
    frame6.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame6, text=frame6.filename)
    my_label.grid()
    return frame6.filename

def muscle():
    os.system('''
    muscle -in {mus} -out {mus}.txt
    '''.format(mus=search2()))


my_btn = Button(frame6, text= "Muscle", command = muscle,padx=40, pady=20,  borderwidth= 5).grid(row=3, column= 0)
##############################################tab 4#############################
frame7 = LabelFrame(tab5, text="Phylogenetic Tree ", padx=150, pady=20)
frame7.pack(padx=20, pady=20)

def search3():
    curr_directory = os.getcwd()
    res5 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your File in newick Formate (.dnd) to draw a phylogenetic tree. make sure that you have to save your photo")
    frame7.filename = filedialog.askopenfilename(initialdir=curr_directory, title="add your sequence")
    my_label = Label(frame7, text=frame7.filename)
    my_label.grid()
    return frame7.filename

def tree():
    my_file = search3()
    from Bio import Phylo
    tree = Phylo.read(my_file, "newick")
    Phylo.draw_ascii(tree)
    tree.rooted = True
    tree.root.color = "salmon"
    tree.clade[0, 1].color = "blue"
    Phylo.draw(tree, branch_labels=lambda c: c.branch_length)

mybtn88= Button(frame7, text="Draw Phylogny", command = tree,  borderwidth= 5).grid()

#################################################tab 6

frame8 = LabelFrame(tab6, text="Sequence Information ", padx=150, pady=20)
frame8.pack(padx=20, pady=20)
from tkinter.ttk import *
from Bio import SeqIO
from Bio.SeqUtils import GC
mylabel = Label(frame8, text= "K-mer Nuclotide Number")
mylabel.grid(row=30, column=1)
combo = Combobox(frame8)

combo['values']= (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

combo.current(2) #set the selected item

combo.grid(column=2, row=30)

def search30():
    res6 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your File in Fasta Formate to get Number of reads, GC content, Number of nuclotides  and headers of your reads.")
    frame8.filename = filedialog.askopenfilename(initialdir="/PycharmProjec/pythonProjects", title="add your sequence")
    my_label = Label(frame8, text=frame8.filename)
    my_label.grid()
    return frame8.filename

def se22():
    mySequence1=search30()
    # count the reads
    myFile = open(f"{mySequence1}.ods", "w")
    count = 0
    for rec in SeqIO.parse(mySequence1, "fasta"):
        count += 1

    myFile.write("Your file contains %i reads" % count)
    myFile.write("\n\n")
    myFile.close()
    ##################################################################################
    # Details of sequences
    # Count letter
    # while loop to counter letter

    myFile = open(f"{mySequence1}.ods", "a")
    myFile.write(f"")
    myFile.write("\n")
    k =combo.get()

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


    for seq_record in SeqIO.parse(mySequence1, "fasta"):
        myFile.write(seq_record.description)
        myFile.write("\n")
        myFile.write(repr(seq_record.seq))
        myFile.write("\n")
        myFile.write("The length of this record is ")
        myFile.write(str(len(seq_record)))  # Length of the sequence
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
        myFile.write("CG content percentage= ")
        myFile.write(str(round(GC(rec.seq), 2)))
        myFile.write("%")
        myFile.write("\n")
        myFile.write("The highest k-mer is => ")
        myFile.write(str(FrequentWords(seq_record.seq, int(k))))
        myFile.write("\n\n\n")

    myFile.write("\n\n")
    myFile.close()
    ################################################################################

mybtn8 = Button(frame8, text= "seq information", command = se22).grid(column= 5, row =30)

#################################################
frame9 = LabelFrame(tab6, text="Reverse Complement ")
frame9.pack(padx=20, pady=20)

def search26():
    res7 = messagebox.showinfo("This is information message", "By clicked okay you wil choose your File in fasta Formate to get reverse complement Sequence")
    frame9.filename = filedialog.askopenfilename(initialdir="/PycharmProjec/pythonProjects", title="add your sequence")
    my_label = Label(frame9, text=frame9.filename)
    my_label.grid()
    return frame9.filename

def Reverse():
    from Bio import SeqIO
    from Bio.Seq import reverse_complement

    mySequence = search26()

    mySeq = open(f"{mySequence}.txt", "w")

    for seq_record in SeqIO.parse(mySequence, "fasta"):
        mySeq.write(seq_record.name)
        mySeq.write("\n")
        mySeq.write(str(reverse_complement(seq_record.seq)))
        mySeq.write("\n\n")
    mySeq.close()
my_btn9 = Button(frame9, text = "Reverse complement sequence", command= Reverse).grid(column= 1, row =2)

window.mainloop()