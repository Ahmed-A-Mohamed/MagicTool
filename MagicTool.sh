# http://genomedata.org/seq-tec-workshop/
echo "Do you need to download the backage? Y or N"
read x
if [ $x == Y ] || [ $x == y ]
then
sudo apt-get update -y
sudo apt-get install python-biopython 
sudo apt-get install python-biopython-doc
sudo apt-get install -y ncbi-entrez-direct 
sudo apt-get install matpotlib 
sudo apt-get install python3-tk 
sudo apt-get install -y clustalw
sudo apt-get install ncbi-blast+
wget https://gist.githubusercontent.com/ozagordi/099bdb796507da8d9426/raw/6ca66616fd545fbb15d94b079e46a7c55edb54c0/blast2sam.py
cd ~
wget https://data.broadinstitute.org/igv/projects/downloads/2.8/IGV_Linux_2.8.13_WithJava.zip
unzip IGV_Linux_2.8.13_WithJava.zip
sudo echo 'export IGV=$HOME/IGV_Linux_2.8.13/igv.sh' >> ~/.bashrc
source ~/.bashrc
sudo apt-get install samtools
fi

echo "Please confirm that your sequence files and script files in the same location"

echo -e "Choose from options (A or B or C ...............)\
\n A. Calling sequence \
\n B. Convert from Fastq to Fasta \
\n C. Pairwise Alignment \
\n D. Multisequence alignment\
\n E. Phylogenetic Tree\
\n F. Reverse Complement\
\n G. information about reads\
\n="
read x

if [ $x == 'A' ] ||  [ $x == 'a' ]
then
  echo "Accession number for ( 1. Gene \ 2. Protein \ 3. protein from gene accession number) = "
  read y
  echo "Attention, if you call many accession numbers, add them with space between every one (ex. NM_001079817.3 NM_000208.4 NM_001322795.2)"
  echo "what is your Accession number?"
  read number
  if [ $y == "1" ]                                                                                        # NM_001079817.3 NM_000208.4
  then
    # for nucleotide
    for i in $number
    do 
    efetch -db nucleotide -id $i -format fasta  >> $i.fasta 
    gedit $i.fasta
    done
  elif [ $y == "2" ]                                                                                     # NP_001355815.1
  then
    # for Protein
    for i in $number
    do
    efetch -db protein -id $i -format fasta >> $i.fasta 
    gedit $i.fasta 
    done
  elif [ $y == "3" ]                                                                                      #  NM_001322795.2
  then
    for i in $number 
    do
    elink -db nuccore  -id $i -target protein | efetch -format fasta > $i.protein.fasta
    gedit $i.protein.fasta 
    done  
    
  else
    echo "PLease, repeat your process"  
  fi  
  echo "Your process is successful"   
###################################################################
elif [ $x == 'B' ] || [ $x == 'b' ] 
then
  echo "Enter your sequence .fastq"
  read seq
  for i in $seq
  do
  cat $seq | paste - - - - | cut -f 1,2 | sed 's/^@/>/' | tr "\t" "\n" >$seq.fasta 
  done

  for filename in $seq.fasta
  do 
  [ -f "$filename" ] || continue
  mv "$filename" "${filename//.fastq/}"
  done
  gedit "${filename//.fastq/}"
  echo "Your process is successful"   

###################################################################
# Pairwise alignment
# 

elif [ $x == 'C' ] || [ $x == 'c' ] 
then
  echo "Pairwise alignment for A.DNA OR B.Protein? Choose A or B"
  read choose
  if [ $choose == 'A' ] || [ $choose == 'a' ]
  then
    echo "Choose from options of pairwise alignment\
     A. your sequence to another sequence\
     B. your sequence to your multisequence database\
     C. your sequence to Human genome 38 database\
     = "
    read z
    if [ $z == 'A' ] || [ $z == 'a' ] 
    then
      echo "Enter first seqences .fasta"
      read first
      echo "Enter second seqences.fasta"
      read second
      blastn -query $first -subject $second -out $first.$second.PairwiseAlignment.txt 
      gedit $first.$second.PairwiseAlignment.txt
    
      echo "PairwiseAlignment.txt and PairwiseAlignment.xml are done, Do you need open alignment by IGV? Y or N"
      read answer
      if [ $answer == "Y" ] || [ $answer == "y" ]
      then
        blastn -query $first -subject $second -outfmt 5 -out $first.$second.PairwiseAlignment.xml
        python3 blast2sam.py $first.$second.PairwiseAlignment.xml > $first.$second.PairwiseAlignment.sam
        samtools view -S -b $first.$second.PairwiseAlignment.sam -o $first.$second.PairwiseAlignment.bam
        samtools sort $first.$second.PairwiseAlignment.bam -o sorted_$first.$second.PairwiseAlignment.bam
        samtools index sorted_$first.$second.PairwiseAlignment.bam
        bash $IGV -g $second sorted_$first.$second.PairwiseAlignment.bam   
      fi 
    
      for filename in $first.$second.PairwiseAlignment.*
      do  
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done
    
      for file in $first.$second.PairwiseAlignment.bam*
      do  
      [ -f "$file" ] || continue
      mv "$file" "${file//.fasta/}"
      done  
      
    ######################################       
    elif [ $z == 'B' ] || [ $z == 'b' ] 
    then
    
      echo "Attention, your database must be indexed, Do you need index your data? Y or N"
      read a
      if [ $a == "Y" ] || [ $a == "y" ]
      then
      makeblastdb -in $data –dbtype nucl
      fi
      
      echo "Confirm that the files (ex. three files) which come from indexing in the same folder with the reference sequences and query sequences" 
      echo "Enter your seqences .fasta"
      read file
      echo "Enter your database .fasta"
      read data
      
      blastn -db $data -query $file -out $data.$file.PairwiseAlignment.txt
      blastn -db $data -query $file -outfmt 6 -out $data.$file.PairwiseAlignmentFilterartion.txt 
      gedit $data.$file.PairwiseAlignment.txt
      gedit $data.$file.PairwiseAlignmentFilterartion.txt
    
      echo "Do you want filteration with e_value?"
      read w
      if [ $w == "Y" ] || [ $w == "y" ]
      then
        echo "Wirte your e_value"
        read e
        awk '{if ($11 <'$e') print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $12}' $data.$file.PairwiseAlignmentFilterartion.txt > Filter.$data.$file.txt
        gedit Filter.$data.$file.txt
      fi  
    
      echo "PairwiseAlignment.txt and PairwiseAlignment.xml are done, Do you need open alignment by IGV? Y or N"
      read answer
      if [ $answer == "Y" ] || [ $answer == "y" ]
      then
        blastn -db $data -query $file -outfmt 5 -out $data.$file.PairwiseAlignment.xml
        python3 blast2sam.py $data.$file.PairwiseAlignment.xml >> $data.$file.PairwiseAlignment.sam
        samtools view -S -b $data.$file.PairwiseAlignment.sam -o $data.$file.PairwiseAlignment.bam
        samtools sort $data.$file.PairwiseAlignment.bam -o sorted_$data.$file.PairwiseAlignment.bam
        samtools index sorted_$data.$file.PairwiseAlignment.bam
        bash $IGV -g $data sorted_$data.$file.PairwiseAlignment.bam
      fi  
    
      for filename in $$data.$file.PairwiseAlignment.*
      do
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done
      for filename in $data.$file.PairwiseAlignment.bam*
      do
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done 
    #########################################
    elif [ $z == 'C' ] || [ $z == 'c' ] 
    then
      echo "Do you need to Download the human genome 38 database (hg38)? Y or N"
      read answer 
      if [ $answer == "Y" ] || [ $answer == "y" ]
      then
         #Download Reference 
         wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
         gunzip hg38.fa.gz
         mv hg38.fa hg38.fasta
         #index your genome
         makeblastdb -in hg38.fasta -dbtype nucl
      fi
      echo "Confirm that the files (three files) which come from indexing in the same folder with the reference sequences and query sequences" 
      echo "Enter your seqences .fasta"
      read file
      blastn -db hg38.fasta -query $file -out hg38.$file.PairwiseAlignment.txt
      blastn -db hg38.fasta -query $file -outfmt 6 -out hg38.$file.PairwiseAlignmentFilterartion.txt
      gedit hg38.$file.PairwiseAlignment.txt
      gedit hg38.$file.PairwiseAlignmentFilterartion.txt
    
      echo "Do you want filteration with e_value?"
      read w
      if [ $w == "Y" ] || [ $w == "y" ]
      then
        echo "Wirte your e_value"
        read e
        awk '{if ($11 <'$e') print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $12}' hg38.$file.PairwiseAlignmentFilterartion.txt > Filter.hg38.$file.txt
        gedit Filter.hg38.$file.txt
      fi
    
      echo "PairwiseAlignment.txt and PairwiseAlignment.xml are done, Do you need open alignment by IGV? Y or N"
      read answer
      if [ $answer == "Y" ] || [ $answer == "y" ]
      then
        blastn -db hg38.fasta -query $file -outfmt 5 -out hg38.$file.PairwiseAlignment.xml
        python3 blast2sam.py hg38.$file.PairwiseAlignment.xml >> hg38.$file.PairwiseAlignment.sam
        samtools view -S -b hg38.$file.PairwiseAlignment.sam -o hg38.$file.PairwiseAlignment.bam
        samtools sort hg38.$file.PairwiseAlignment.bam -o sorted_hg38.$file.PairwiseAlignment.bam
        samtools index sorted_hg38.$file.PairwiseAlignment.bam
        bash $IGV -g hg38.fasta sorted_hg38.$file.PairwiseAlignment.bam

      fi   
    
      for filename in hg38.$file.PairwiseAlignment.*
      do  
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done
    
      for filename in hg38.$file.PairwiseAlignment.bam*
      do  
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done

    else
      echo 'Please repeat the process'
    fi

  echo "Your process is successful"
#######################################################################################
#######################################################################################  
  elif [ $choose == 'B' ] || [ $choose == 'b' ]
  then
    echo "Choose from options of pairwise alignment\
     A. your sequence to another sequence\
     B. your sequence to your database\
     C. your sequence to Human protein database\
     = "
    read z
    if [ $z == 'A' ] || [ $z == 'a' ] 
    then
      echo "Enter first seqences .fasta"
      read first
      echo "Enter second seqences.fasta"
      read second
      blastp -query $first -subject $second -out $first.$second.PairwiseAlignment.txt 
      gedit $first.$second.PairwiseAlignment.txt 

      for filename in $first.$second.PairwiseAlignment.*
      do  
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done
    ######################################       
    elif [ $z == 'B' ] || [ $z == 'b' ] 
    then
      echo "Attention, your database must be indexed, Do you need index your data? Y or N"
      read a
      if [ $a == "Y" ] || [ $a == "y" ]
      then
      makeblastdb -in $data –dbtype prot
      fi
    
      echo "Confirm that the files (ex. three files) which come from indexing in the same folder with the reference sequences and query sequences" 
      echo "Enter your seqences .fasta"
      read file
      echo "Enter your database .fasta"
      read data
    
      blastp -db $data -query $file -out $data.$file.PairwiseAlignment.txt
      blastp -db $data -query $file -outfmt 6 -out $data.$file.PairwiseAlignmentFilterartion.txt 
      gedit $data.$file.PairwiseAlignment.txt
      gedit $data.$file.PairwiseAlignmentFilterartion.txt
    
      echo "Do you want filteration with e_value?"
      read w
      if [ $w == "Y" ] || [ $w == "y" ]
      then
        echo "Wirte your e_value"
        read e
        awk ' {if ($11<'$e') print $1 "\t" $2 "\t" $11}' $data.$file.PairwiseAlignmentFilterartion.txt > filter.$data.$file.txt
        gedit filter.$data.$file.txt
      fi  
    
    
      for filename in $$data.$file.PairwiseAlignment.*
      do
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done

    #########################################
    elif [ $z == 'C' ] || [ $z == 'c' ] 
    then
      echo "Do you need to Download the human protein database ? Y or N"
      read answer 
      if [ $answer == "Y" ] || [ $answer == "y" ]
      then
         #Download Reference 
         wget ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/mRNA_Prot/human.1.protein.faa.gz
         gunzip human.1.protein.faa.gz
         mv human.1.protein.faa human.protein.fasta
         #index your protein
         makeblastdb -in human.protein.fasta -dbtype prot
         
      fi
      echo "Confirm that the files (ex. three files) which come from indexing in the same folder with the reference sequences and query sequences" 
      echo "Enter your seqences .fasta"
      read file
      blastp -db human.protein.fasta -query $file -out human.protein.$file.PairwiseAlignment.txt
      blastp -db human.protein.fasta -query $file -outfmt 6 -out human.protein.$file.PairwiseAlignmentFilterartion.txt
      gedit human.protein.$file.PairwiseAlignment.txt
      gedit human.protein.$file.PairwiseAlignmentFilterartion.txt
    
      echo "Do you want filteration with e_value?"
      read w
      if [ $w == "Y" ] || [ $w == "y" ]
      then
        echo "Wirte your e_value"
        read e
        awk ' {if ($11<'$e') print $1 "\t" $2 "\t" $11}' human.protein.$file.PairwiseAlignmentFilterartion.txt > filter.human.protein.$file.txt 
        gedit filter.human.protein.$file.txt
      fi  
    
      for filename in human.protein.$file.PairwiseAlignment.*
      do  
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done

    else
      echo 'Please repeat the process'
    fi

  echo "Your process is successful"
     
  else
    echo 'Please repeat the process'
  fi    
###################################################################
elif [ $x == 'D' ] || [ $x == 'd' ] 
then
  echo "Do you need (A. ClustelW) OR (B. Muscle)?"
  read answer
  if [ $answer == "A" ] || [ $answer == "a" ]
  then  
    echo "Enter your sequence in .fasta format!"
    read seq
    export seq
    python3 MultisequenceAlignment.py
    gedit *.aln 
    gedit *.dnd
  
  elif [ $answer == "B" ] || [ $answer == "b" ]
  then 
    echo "Enter your sequence .fasta format"
    read seq
    muscle -in $seq -out $seq.txt
    
    for filename in $seq.txt
      do  
      [ -f "$filename" ] || continue
      mv "$filename" "${filename//.fasta/}"
      done
      gedit "${filename//.fasta/}"
   
   else
     echo "Please. repeat your process"   
   fi  

echo "Your process is successful"
###################################################################
elif [ $x == 'E' ] || [ $x == 'e' ] 
then
  python3 Phylogenetic.py
  
echo "Your process is successful"  
###################################################################
elif [ $x == 'F' ] || [ $x == 'f' ] 
then
  echo "Enter your sequence in .fasta format!"
  read seq
  export seq
  python3 'Reverse complement.py'
  gedit reverse_complement.$seq.txt

echo "Your process is successful" 
###################################################################
elif [ $x == 'G' ] || [ $x == 'g' ] 
then

  python3 Infotmation_sequence.py
  
echo "Your process is successful" 
###################################################################

else
  echo 'Please repeat the process'
fi
echo "Do you need another process? Y or N"
read another
if [ $another == "Y" ] || [ $another == "y" ]
then 
  bash MagicTool.sh
else
  echo "Thanks for using this tool"
fi 

