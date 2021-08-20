# Exaiptasia diaphana transcriptomics analysis workflow


### Context: 
The data used is from the paper “Seneca et al. 2020 BMC Genomics “Gene expression kinetics of Exaiptasia  pallida innate immune response to Vibrio parahaemolyticus infection”

BIOPROJECT: PRJNA665134


I have downloaded the fastq files for the RNA Seq data from the T= 12 hr time points from the triplicate samples of what I assume is: the treatment with Vibrio Para. And the control (Accession numbers: SRR12696075 SRR12696076 SRR12696077 SRR12696079 SRR12696080 SRR12696081)


I will be using this small data set to try to run some basic transcriptomic analyses and see what comes out of this mess. 


- GTF - Gene Annotation for E. diaphana downloaded from NCBI Genome 
- Genome of E. diaphana downloaded from NCBI SRA (Baumgarten et al. 2015 PNAS)




## workflow

Tutorial: 
https://biocorecrg.github.io/RNAseq_course_2019/multiqc.html


Directories to be generated 
1. data  
   1a. qc #where qc data will be stored  
   1b. annotation #where annotations and reference genome is stored  
   1c. indices #where generated index genome is stored  
   1d. alignments #where generated alignments are stored  
   1e. deseq2 


_Notes_  
`zcat` to inspect/open zipped filed  
`cat` used to inspect/ open unzipped filed


### 1. Download genome and annotations 
Keep the genome and annotation as .gz files and only unzip for temporary use. Remove the unzipped files once done with them. They are too large for long term storage 

.fa and .fasta are the same file types, can be renamed to work with specific programs.  

### 2. Inspect genome and annotations 

```markdown
mkdir annotations 
```

#where GTF file and genome.fa can live





### 3. FASTQ Transcript sequence files

1. place FASTQ files in data folder

2. rename fastq files with lookup table 

seeking to replace accession numbers with sample id’s

loop with lookup table? ( i could do this in r with a data frame easily, but dont know how to do in BASH/ Python

```markdown
for file in $(ls *.fastq	| sed -r 's/
```

----- tony's idea --------

As long as the files are in the same order as the fastq files then you could do something like below.
Make a directory for the forward, then move them there, same for reverse.
Then get just column two of the sra_lookup.txt file (sra_lookup_two.txt).

then something like this:

for file in *.fastq; do read line; mv -v "${file}" "${line}"; done < sra_lookup_two.txt

Then repeat for the reverse. 

-------------------------



### 4. fastQC and fastq screen 
I have downloaded the interactive application, but haven’t run the command line version. Need this to use multiqc report compilation later

-------------------------


### 5. STAR - sequence alignment 

*apparently requires a great deal of RAM to run, my computer seems to work though*

Compiled STAR unix executable file being stored locally `users/local_computation` (not in my usual cloud storage)

#### 1. Prep directories for analysis


```markdown
cd local_computation/STAR-2.7.9a/bin 
#set wd for the STAR executable files

export PATH=$PATH:$PWD 
# add executable binaries to PATH variable (doesn’t require absolute paths) 

echo $PATH 
# to check appending successful 

chmod +x ./STAR 
#when first running code, needed to give permission to access the STAR file 

cd ~/Library/Mobile\ Documents/com~apple~CloudDocs/WorkDocuments/projects/exaiptasia_rnaseq_practice/data 
#change WD back to directory with the reference genome MAY REQUIRE UPDATES IF PATH CHANGES
```


#### 2. Make index - pretty quick step

it might be best to put all this into its own directory to keep files separate…

mkdir indices/ExdiChr 
#where indexed exaiptasia diaphana (Exdi) chromosome can live  


version 1.1
```markdown
STAR --runMode genomeGenerate \
 --runThreadN 8 \
 --genomeDir indices/ExdiChr \
 --genomeFastaFiles annotations/Aiptasia_genome_1-1_genomic.fa \
 --genomeSAindexNbases 6 \
 --sjdbGTFfile annotations/exaiptasiaDiaphana_genomic.gtf \
 --sjdbOverhang 98 \
--outFileNamePrefix exdi_index

#Code ANNOTATION
#STAR  --runMode genomeGenerate \
# --runThreadN 8 \ #how many cores are available for processing, more is faster
# --genomeDir ./ \ # this should map to where the genome is stored. If WD is already this location, `./ ` indicates look in WD
# --genomeFastaFiles Aiptasia_genome_1-1_genomic.fa \  #should map to reference genome 
# --genomeSAindexNbases 6 \ #this is used because the genome I am working with is very small, number can vary
# --sjdbGTFfile exaiptasiaDiaphana_genomic.gtf \ #merge with available gene annotations 
# --sjdbOverhang 98 \  #read length -1 for something to do with splicing 
```


#### 3. align RNASEQ reads to indexed genome - much more time intensive step

Check file capacity 
```markdoen
ulimit -n 
#check number of files allowable to be open at once, if lower than 10 000, increase tp 10 000

ulimit -n 10000 
#increase number of files open to 10 000
```

create directory for aligned transcripts
```markdown
mkdir alignments 
# make a directory to store generated alignment data - check paths in use files to make sure its being used 
```




Version 1.1 I think STAR needs to be run iteratively with each sample type to create an output file for each one, not run with a big string of fasta files after it. This code is an attempt to use a loop to run though a subset of fastq files 

** there are options to load and retain indexed genome, cant figure out how to do it yet **

``` markdown
STAR --genomeLoad LoadAndExit --genomeDir index 
#cant do this until shared memory increased, without this command, the loop for alignment will be slowed considerably if using a larger genome. Works ok for exaiptasia

for i in $(ls *.fastq | sed -r 's/_[12][.]fastq//' | uniq)
do 
STAR --genomeDir indices/ExdiChr \
--readFilesIn ${i}_1.fastq ${i}_2.fastq \
--runThreadN 8 \
--outFileNamePrefix alignments/Exdi_STAR_$i. \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts ; done


#this command shows what is happening when the loop is run, useful for troubleshooting **

# for i in $(ls *.fastq | sed -r 's/_[12][.]fastq//' | uniq)
# do 
# echo STAR --genomeDir indices/ExdiChr \
# --readFilesIn ${i}_1.fastq ${i}_2.fastq \
# --runThreadN 8 \
# --outFileNamePrefix alignments/Exdi_STAR_$i. \
# --outSAMtype BAM SortedByCoordinate \
# --quantMode GeneCounts ; done
```

-------------------------



#### 4. Inspect reads per gene file - figure out which strandedness in transcripts
find out strandedness of RNA sequencing preparation, see [Biocore RNAseq course] https://biocorecrg.github.io/RNAseq_course_2019/

```markdown
more alignments/Exdi_STAR_T-12h-1.ReadsPerGene.out.tab
#display reads per gene in LONG list


grep -v "N_" alignments/Exdi_STAR_T-12h-1.ReadsPerGene.out.tab | awk '{unst+=$2;forw+=$3;rev+=$4}END{print unst,forw,rev}'
#count reads per strand/ gene for one .tab file 
```

-------------------------



### 6. samtools - .bam => .sam => .cram
```markdown
cd local_computation/samtools-1.13 
#change directory to unix executable directory


export PATH=$PATH:$PWD  
#add. Samtools to path

cd ~/Library/Mobile\ Documents/com~apple~CloudDocs/WorkDocuments/projects/exaiptasia_rnaseq_practice/data 
#change WD back to directory with the reference genome MAY REQUIRE UPDATES IF PATH CHANGES
```



#### 1. inspect .bam files

```markdown
for i in $(ls alignments/Exdi_STAR_*.Aligned.sortedByCoord.out.bam | uniq)
do 
samtools view -h $i | head -n 10 ; done

#this shows the ten first rows of the .bam files in the directory


#echo check
#for i in $(ls alignments/Exdi_STAR_*.Aligned.sortedByCoord.out.bam | uniq)
#do 
#echo samtools view -h $i | head -n 10 ; done
```

#### 2. .BAM => .SAM

```markdown
for i in $(ls alignments/*.bam | sed -r 's/[.]bam//' | uniq)
do 
samtools view -h ${i}.bam > ${i}.sam ; done


#echo check
#for i in $(ls alignments/*.bam | sed -r 's/[.]bam//' | uniq)
#do 
#echo samtools view -h ${i}.bam > ${i}.sam ; done
```


#### 3. .SAM => .CRAM
```markdoen
samtools faidx annotations/Aiptasia_genome_1-1_genomic.fa

for i in $(ls alignments/*.sam | sed -r 's/[.]sam//' | uniq)
do
samtools view -C ${i}.sam -T annotations/Aiptasia_genome_1-1_genomic.fa > ${i}.cram ; done
```

-------------------------


### 7. multiQC - check quality scores (havent done this yet)  
Generate report for mapping quality, direct to the directory generated for fastqc reports 


#conda install -c bioconda -c conda-forge multiqc #install multiqc
#Set to directory with quality control data 

cd qc

multiqc . #run multiqc analysis and generate html report 

-------------------------


### 8. DESeq2 - DGE (differential gene expression)

http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html
https://hbctraining.github.io/DGE_workshop/lessons/04_DGE_DESeq2_analysis.html

will need to format gene data so that it is compatible with deseq (matrix: rows as genes, columns as samples)



#### 1. Make matrix with gene expression data for each sample 

```markdown
mkdir deseq2 
#make directory for deseq2 analysis.  
```

 **identify strandedness of RNA seq analysis - see tutorial (https://biocorecrg.github.io/RNAseq_course_2019/differential_expression.html). Whatever columns of the counts represent the correct reads, change below code to extract correct column **


```markdown
paste alignments/*ReadsPerGene.out.tab | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > alignments/tmp

# this tmp file adds an extra column to the right of the data, could do a second bit of code to only select the relevant columns that contain text...
# retrieve the 4th column of each "ReadsPerGene.out.tab" file + the first column that contains the gene IDs
```
```markdown
# this code doesn't work with my test data because each gene ID has an underscore (AC249_AIPGENE13054) which is then removed by the grep -v "_" command and then it has only 3 columns which prevents the awk script from pulling the $4th column and populating the tmp file. If the gene ID was formatted like the example data "ENSG00000223972.5", there would be no problem. Keeping the example script incase the need artises. In my situation, I just need to delete the rows with the "N_unmapped, N_multimapping, N_noFeature and N_ambigious" before trying to ID genes present here.
# paste *ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp
```

#### 1a. remove first four rows of unwanted text  
manually remove the first four lines of the tmp file because I cant be bothered to com up with a script to do it RN. Maybe with `awk`... or just refine the previous awk script that made the tmp file


#### 1b. add column headers 

```markdown
ls alignments/*ReadsPerGene.out.tab | awk 'BEGIN{ORS="";print "gene name\t"}{print $0"\t"}END{print "\n"}'| sed 's/[.]ReadsPerGene.out.tab//g' > deseq2/raw_counts_Exdi_matrix.txt ; cat alignments/tmp >> deseq2/raw_counts_Exdi_matrix.txt

# add header: "gene_name" + the name of each of the counts file
```


#### 2. Make annotation file a .csv - problem step 

** this file will require a bit of new code to make the proper file. Will compare between full_data file and my exaiptasia file to decide what details should be included in the tx2gene file. **

`cat annotations/exaiptasiaDiaphana_genomic.gtf | awk -F "\t" 'BEGIN{OFS="\t"}{if($3=="transcript"){split($9, a, "\""); print a[4],a[2],a[8]}}' > deseq2/tx2gene.exdi.csv`

#### 3. Make sample sheet - code not run

done manually in excel (be careful with forward slashes, R will turn headers with them into periods... incentive to remove `annotations/` from the begining of all my sample headers...
`cat <(echo -e "SampleName\tFileName\tTime\tDexamethasone") <(paste <(ls counts_4thcol | cut -d"_" -f1-3) <(ls counts_4thcol) <(ls counts_4thcol | cut -d"_" -f2 | awk '{print "t"$0}') <(printf '100nM\n%.0s' {1..6})) > sample_sheet_Exdi.txt`



#### 4. Open Deseq2 in R 

```markdown
R #open R

getwd

setwd("~/deseq2")
```



















