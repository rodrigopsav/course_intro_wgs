# Tools for whole genomic sequencing analysis

## About the course

The main goals of this workshop are:

**(1)** to present the most used file formats to keep the WGS information.
**(2)** to give an overview of each workflow step to get the variants from the raw fastq files.
**(3)** to use bioinformatic tools to analyse sequencing data.
**(4)** to present IVDP: Integrated Variant Discovery Pipeline

[MSU HPCC](#msu-hpcc)   
[Download files](#download-files)   
[Download and install IVDP](#download-and-install-ivdp)   
[Install IVDP dependencies](#install-ivdp-dependencies)   


[Practice 1: file formats](#practice-1-file-formats)   
* [Searching for HPCC modules](#searching-for-HPCC-modules)
* [Loading HPCC modules](#loading-HPCC-modules)
* [FASTA format](#fasta-format)
* [Indexing FASTA with samtools](#indexing-fasta-with-samtools)
* [FASTQ format](#fastq-format)
* [SAM / BAM formats](#sam--bam-formats)
* [VCF format](#vcf-format)
* [GTF format](#gtf-format)


[Practice 2: WGS workflow](#practice-2-wgs-workflow)   
* [Creating directories](#creating-directories)
* [Indexing reference genome](#indexing-reference-genome)
* [Selecting sample](#selecting-sample)
* [QC of FASTQ files](#qc-of-fastq-files)
* [Trimming reads](#trimming-reads)
* [QC of trimmed fastq files](#qc-of-trimmed-fastq-files)
* [Alignment](#alignment)
* [Mark duplicated reads](#mark-duplicated-reads)
* [Base quality score recalibration (BQSR)](#base-quality-score-recalibration-bqsr)
* [VCF calling](#vcf-calling)


[Practice 3: IVDP](#practice-3-ivdp)   
* [Running IVDP](#running-ivdp)


## MSU HPCC

### Login on HPCC
```
ssh -XY user@hpcc.msu.edu
password:
ssh dev-intel18
```

### HPCC tutorials
[hpcc wiki](https://wiki.hpcc.msu.edu/display/ITH/High+Performance+Computing+at+ICER)   
[open OnDemand](https://wiki.hpcc.msu.edu/display/ITH/Open+OnDemand)   
[login OnDemand](https://ondemand.hpcc.msu.edu/)   

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


## Download files
```
cd
git clone https://github.com/rodrigopsav/course_intro_wgs.git
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


## Download and install IVDP
```
# Download IVDP
cd ~/
mkdir softwares
cd softwares
git clone https://github.com/rodrigopsav/IVDP.git
wait


# Install IVDP
cd IVDP/install_ivdp_dependencies
./install_ivdp_dependencies.sh -d ~/softwares
cd ..
wait
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


## Practice 1: file formats


### Searching for HPCC modules
```
module spider samtools
module spider SAMtools/1.10
```

### Loading HPCC modules
```
module purge
module load GCC/8.3.0
module load HTSlib/1.10.2
module load SAMtools/1.10
module load BCFtools/1.10
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### FASTA format
```
cd course_intro_wgs/files/examples

# Print fasta file page by page
cat genome.fa | less

# Search strings with ">"
cat genome.fa | grep ">"

# Count number of strings with ">"
cat genome.fa | grep ">" | wc -l

# Search string with number 2
cat genome.fa | grep ">2"

# Search string with number 2 and get the first result 
cat genome.fa | grep -m 1 ">2"

# Print 3 lines after the result of "grep" search
cat genome.fa | grep -m 1 -A 3 ">2"

# Print 3 lines after and 7 lines before the result of "grep" search
cat genome.fa | grep -m 1 -A 3 -B 7 ">2"
```
-m 1, --max-count, stop reading a file after 1 matching line   
-A, show lines after grep   
-B, show lines before grep   

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Indexing FASTA with samtools
```
samtools faidx genome.fa
cat genome.fa.fai
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### FASTQ format
```
zcat example9_S2_R1_001.R1.fastq.gz | less
```

[Naming convention rules](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm) for FASTQ files involves:

▶ SampleName — The sample name provided in the sample sheet. If a sample name is not provided, the file name includes the sample ID, which is a required field in the sample sheet and must be unique.   
▶ S1 — The sample number based on the order that samples are listed in the sample sheet starting with 1. In this example, S1 indicates that this sample is the first sample listed in the sample sheet.   
▶ L001—The lane number.   
▶ R1—The read. In this example, R1 means Read 1. For a paired-end run, there is at least one file with R2 in the file name for Read 2. When generated, index reads are I1 or I2.   
▶ 001—The last segment is always 001.   

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### SAM / BAM formats
```
# View sam file
samtools view example.sam | less

# Count number of reads
samtools view example.sam | wc -l

# -h prints sam file with header
samtools view -h example.sam | less

# -H prints only the header
samtools view -H example.sam

# Sort sam and convert to bam file
samtools sort example.sam -o example_sorted.bam

# Index bam file
samtools index example_sorted.bam

# Print bam file
samtools view example_sorted.bam | less

# Print reads from chromosome 1
samtools view example_sorted.bam 1 | wc -l

# Print reads from chromosome 2
samtools view example_sorted.bam 2 | wc -l

# Count number of reads by region
samtools view example_sorted.bam 1:1-100 | wc -l

# Print reads by region
samtools view example_sorted.bam 1:1-100 | awk '{print $1, $4, $5, $6, $10}'
samtools view example_sorted.bam 1:1-100 | awk '{print $1, $4, $5, $6, length($10)}'

# Bam file stats
samtools stats example_sorted.bam > exampleBam.stats
less exampleBam.stats
cat exampleBam.stats | grep "^SN"

# More stats
samtools flagstat example_sorted.bam

# Coverage depth
samtools coverage example_sorted.bam | less
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### VCF format
```
# Compress vcf file
bgzip -@ 16 example.vcf

# Index vcf file
bcftools index -f -t --threads 16 example.vcf.gz

# Print vcf file
bcftools view example.vcf.gz | less

# -h prints only vcf's header
bcftools view -h example.vcf.gz

# -H prints vcf without header
bcftools view -H example.vcf.gz | less

# Query information: sample list
bcftools query -l example.vcf.gz

# Query information: map
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' example.vcf.gz | less

# Query information: genotypes
bcftools query -f '%CHROM\t%POS\t%ID[\t%GT]\n' example.vcf.gz | less

# Query information: coverage depth by genotype
bcftools query -f '%CHROM\t%POS\t%ID[\t%DP]\n' example.vcf.gz | less

# Query information: maf
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/MAF\n' example.vcf.gz | less
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### GTF format
```
# Download gtf
wget http://ftp.ensembl.org/pub/release-104/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.104.gtf.gz
gunzip Bos_taurus.ARS-UCD1.2.104.gtf

# Printing gtf
less Bos_taurus.ARS-UCD1.2.104.gtf
awk -F"\t" '{print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf  | less

# Printing without header
awk -F"\t" '{print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf  | grep -v '^#' | less

# Searching with awk
awk -F"\t" '$3 == "gene" {print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf  | grep -v '^#' | less
awk -F"\t" '$3 == "gene" && $1 == 3 {print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf  | grep -v '^#' | less
awk -F"\t" '$3 == "gene" && $1 == 3 && $4 -gt 100000 && $5 -lt 500000 {print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf  | grep -v '^#' | less
awk -F"\t" '$3 == "gene" && $1 == 3 && $4 >= 30000000 && $5 <= 32000000 {print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf | less
awk -F"\t" '$3 == "gene" && $1 == 3 && $4 >= 30000000 && $5 <= 32000000 {print $1, $3, $4, $5}' Bos_taurus.ARS-UCD1.2.104.gtf | wc -l
awk -F"\t" '$3 == "gene" && $1 == 3 && $4 >= 30000000 && $5 <= 32000000 {print $0}' Bos_taurus.ARS-UCD1.2.104.gtf | less

# Searching with grep
grep -v "^#" Bos_taurus.ARS-UCD1.2.104.gtf | grep -i "dgat" | awk '$3 == "gene" {print}'
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


## Practice 2: WGS workflow

### Creating directories
```
# Creating directories
cd
cd course_intro_wgs/files

mkdir -p output/01_qc
mkdir -p output/02_trim
mkdir -p output/03_qcTrimmed
mkdir -p output/04_align
mkdir -p output/05_mdup
mkdir -p output/06_bqsr
mkdir -p output/07_call
mkdir -p output/tmp
mkdir -p output/index
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Indexing reference genome

```
# Loading modules
module purge
module load ifort/2018.3.222-GCC-7.3.0-2.30  impi/2018.3.222
module load bwa-mem2/2.0pre2
module load SAMtools/1.9
module load BCFtools/1.9
module load picard/2.22.1-Java-11

# Index fasta: creating .fai
samtools faidx refGen/reference.fa

# Index vcf: creating .tbi
bcftools index -t refGen/variants.vcf.gz

# Create genome dictionary
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
   R=refGen/reference.fa \
   O=refGen/reference.dict

# Index: for bwa-mem2
bwa-mem2 index -p output/index/refGenIdx refGen/reference.fa
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Selecting sample

```
SAMPLENAME=sample1
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### QC of FASTQ files

```
module load FastQC/0.11.7-Java-1.8.0_162

fastqc -o output/01_qc \
   -t 16 \
   -dir output/tmp \
data/"$SAMPLENAME"_R1.fastq.gz data/"$SAMPLENAME"_R2.fastq.gz
 ```
   
<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Trimming reads

```
module load Trimmomatic/0.39-Java-11

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 -phred33 \
   -summary output/02_trim/"$SAMPLENAME".summary \
   data/"$SAMPLENAME"_R1.fastq.gz data/"$SAMPLENAME"_R2.fastq.gz \
   output/02_trim/"$SAMPLENAME"_R1_paired.fastq.gz output/02_trim/"$SAMPLENAME"_R1_unpaired.fastq.gz \
   output/02_trim/"$SAMPLENAME"_R2_paired.fastq.gz output/02_trim/"$SAMPLENAME"_R2_unpaired.fastq.gz \
   ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:8:true \
   LEADING:15 TRAILING:15 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:35 CROP:250
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### QC of trimmed fastq files

```
fastqc -o output/03_qcTrimmed \
   -t 16 \
   -dir output/tmp \
   output/02_trim/"$SAMPLENAME"_*_paired.fastq.gz
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Alignment

```
# Loading modules
module load bwa-mem2/2.0pre2
module load SAMtools/1.9

# Alignment
bwa-mem2 mem -t 16 \
   -o output/04_align/"$SAMPLENAME"_trim_bwa2.sam \
   output/index/refGenIdx \
   output/02_trim/"$SAMPLENAME"_R1_paired.fastq.gz output/02_trim/"$SAMPLENAME"_R2_paired.fastq.gz

# Convert sam to bam and sort
samtools sort output/04_align/"$SAMPLENAME"_trim_bwa2.sam -o output/04_align/"$SAMPLENAME"_trim_bwa2_sorted.bam
samtools index output/04_align/"$SAMPLENAME"_trim_bwa2_sorted.bam
rm output/04_align/"$SAMPLENAME"_trim_bwa2.sam

# Add read group
module load picard/2.22.1-Java-11

java -Xmx20g -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
   I=output/04_align/"$SAMPLENAME"_trim_bwa2_sorted.bam \
   O=output/04_align/"$SAMPLENAME"_trim_bwa2_sorted_RG.bam \
   RGID="$(zcat -f data/"$SAMPLENAME"_R1.fastq.gz | head -1 | awk -F ":" '{ print $3}')" \
   RGLB=lib1 \
   RGPL=illumina \
   RGPU="$(zcat -f data/"$SAMPLENAME"_R1.fastq.gz | head -1 | awk -F ":" '{ print $3}')" \
   RGSM="$SAMPLENAME" \
   VALIDATION_STRINGENCY=LENIENT
samtools index output/04_align/"$SAMPLENAME"_trim_bwa2_sorted_RG.bam
rm output/04_align/"$SAMPLENAME"_trim_bwa2_sorted.bam*
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Mark duplicated reads

```
module load Sambamba/0.7.1

sambamba markdup -t 8 --overflow-list-size 1000000 \
   --hash-table-size 1000000 --tmpdir=output/index \
   output/04_align/"$SAMPLENAME"_trim_bwa2_sorted_RG.bam \
   output/05_mdup/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup.bam
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Base quality score recalibration (BQSR)

```
module purge
module load GNU/6.4.0-2.28  OpenMPI/2.1.2-CUDA
module load GATK/4.1.4.1-Python-3.6.4
module load R/3.5.0-X11-20180131


gatk --java-options '-Djava.io.tmpdir=output/tmp -Xmx20g' BaseRecalibrator \
   -R refGen/reference.fa \
   -I output/05_mdup/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup.bam \
   -known-sites refGen/variants.vcf.gz \
   -O output/06_bqsr/"$SAMPLENAME"_recalibration_run1.table \
   --tmp-dir output/tmp

# Apply BQSR
gatk --java-options '-Djava.io.tmpdir=output/tmp -Xmx20g' ApplyBQSR \
   -R refGen/reference.fa \
   -I output/05_mdup/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup.bam \
   --bqsr-recal-file output/06_bqsr/"$SAMPLENAME"_recalibration_run1.table \
   -O output/06_bqsr/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup_bqsr.bam \
   --tmp-dir output/tmp

# (Optional) BaseRecalibration: run2
gatk --java-options '-Djava.io.tmpdir=output/tmp -Xmx20g' BaseRecalibrator \
   -R refGen/reference.fa \
   -I output/06_bqsr/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup_bqsr.bam \
   -known-sites refGen/variants.vcf.gz \
   -O output/06_bqsr/"$SAMPLENAME"_recalibration_run2.table \
   --tmp-dir output/tmp

# (Optional) Plot a single recalibration table
gatk AnalyzeCovariates \
   -bqsr output/06_bqsr/"$SAMPLENAME"_recalibration_run1.table \
   -plots output/06_bqsr/"$SAMPLENAME"_AnalyzeCovariates_run1.pdf

# (Optional) Plot "before" (first pass) and "after" (second pass) recalibration tables to compare them
gatk AnalyzeCovariates \
   -before output/06_bqsr/"$SAMPLENAME"_recalibration_run1.table \
   -after output/06_bqsr/"$SAMPLENAME"_recalibration_run2.table \
   -plots output/06_bqsr/"$SAMPLENAME"_AnalyzeCovariates_compare.pdf
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### Run sample2
```
SAMPLENAME=sample2
```
Run from [QC of FASTQ files](#qc-of-fastq-files) to Base quality score recalibration (BQSR)

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


### VCF calling

```
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6
module load BCFtools/1.9


### HaplotypeCaller with GVCF mode
for SAMPLENAME in sample1 sample2; do
   gatk --java-options '-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=4' HaplotypeCaller \
       -R refGen/reference.fa \
       -I output/06_bqsr/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup_bqsr.bam  \
       -pairHMM LOGLESS_CACHING \
       --native-pair-hmm-threads 8 \
       -O output/07_call/"$SAMPLENAME".g.vcf \
       -ERC GVCF \
       --sample-ploidy 2 \
       --tmp-dir output/tmp &
done


### GenomicsDBImport
# Creating sample list (gvcf files)
cd output/07_call/
ls *g.vcf > input_gvcf.list

rm samples.txt
for input in $(cat input_gvcf.list); do
   bcftools query -l $input >> samples.txt
done
wait

paste -d"\t" samples.txt input_gvcf.list > cohort.sample_map
wait


# Combining GVCFs
rm -r DB_*
for CHROM in 1 2; do
   gatk --java-options '-Xmx16g' GenomicsDBImport \
   --genomicsdb-workspace-path DB_"$CHROM" \
   --sample-name-map cohort.sample_map \
   --batch-size 50 \
   -L $CHROM \
   --reader-threads 8 \
   --tmp-dir ../tmp &
done

### GenotypeGVCFs: convert gvcf to vcf
for CHROM in 1 2; do
   gatk --java-options '-Xmx16g' GenotypeGVCFs \
   -R ../../refGen/reference.fa \
   -V gendb://DB_"$CHROM" \
   -G StandardAnnotation \
   -O bwa2_gatk4_"$CHROM"_raw.vcf \
   -L $CHROM \
   --sample-ploidy 2 \
   --tmp-dir ../tmp &
done


### Merging vcfs by chromosome
# Extracting head of vcf
cat $(ls bwa2_gatk4_*_raw.vcf | head -1) | grep "^#" > bwa2_gatk4_raw.vcf
wait

# Joinning vcfs
for CHROM in 1 2; do
   grep -v '^#' bwa2_gatk4_"$CHROM"_raw.vcf >> bwa2_gatk4_raw.vcf
done
wait

bgzip -f -@ 16 bwa2_gatk4_raw.vcf
bcftools index -f -t --threads 16 bwa2_gatk4_raw.vcf.gz
wait

rm -r input_gvcf.list samples.txt cohort.sample_map DB_* bwa2_gatk4_*_raw.*
cd ../..
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>


## Practice 3: IVDP

### Running IVDP

```
### Run on a local server
# Example 2
./ivdp.sh -p params/pe_ex2.txt

# Example 3b
./ivdp.sh -p params/pe_ex3b.txt

# Example 10
./ivdp.sh -p params/se_ex10.txt

# Example 12
./ivdp.sh -p params/se_ex12.txt


### Run on HPCC with slurm
# Example 2
./ivdp.sh -p params/pe_ex2.txt -c params/configSlurm.txt

# Example 3b
./ivdp.sh -p params/pe_ex3b.txt -c params/configSlurm.txt

# Example 10
./ivdp.sh -p params/se_ex10.txt -c params/configSlurm.txt

# Example 12
./ivdp.sh -p params/se_ex12.txt -c params/configSlurm.txt

```
**Note:** To open html and pdf files with command line, use:
```
evince file.pdf
firefox file.html
```

<div align="right">
    <b><a href="#about-the-course">↥ back to top</a></b>
</div>

