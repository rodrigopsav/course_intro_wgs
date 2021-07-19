#-------------------------------------#
#     Practice 2: WGS workflow        #
#-------------------------------------#


### Creating directories ####
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

  
### Indexing reference genome ####
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


#### Selecting sample ####
SAMPLENAME=sample1

  
### QC of FASTQ files ####
module load FastQC/0.11.7-Java-1.8.0_162

fastqc -o output/01_qc \
-t 16 \
-dir output/tmp \
data/"$SAMPLENAME"_R1.fastq.gz data/"$SAMPLENAME"_R2.fastq.gz


### Trimming reads ####
module load Trimmomatic/0.39-Java-11

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 -phred33 \
-summary output/02_trim/"$SAMPLENAME".summary \
data/"$SAMPLENAME"_R1.fastq.gz data/"$SAMPLENAME"_R2.fastq.gz \
output/02_trim/"$SAMPLENAME"_R1_paired.fastq.gz output/02_trim/"$SAMPLENAME"_R1_unpaired.fastq.gz \
output/02_trim/"$SAMPLENAME"_R2_paired.fastq.gz output/02_trim/"$SAMPLENAME"_R2_unpaired.fastq.gz \
ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:8:true \
LEADING:15 TRAILING:15 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:35 CROP:250

  
### QC of trimmed fastq files ####
fastqc -o output/03_qcTrimmed \
-t 16 \
-dir output/tmp \
output/02_trim/"$SAMPLENAME"_*_paired.fastq.gz


### Alignment ####
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


### Mark duplicated reads ####
module load Sambamba/0.7.1

sambamba markdup -t 8 --overflow-list-size 1000000 \
--hash-table-size 1000000 --tmpdir=output/index \
output/04_align/"$SAMPLENAME"_trim_bwa2_sorted_RG.bam \
output/05_mdup/"$SAMPLENAME"_trim_bwa2_sorted_RG_mdup.bam

  
### Base quality score recalibration (BQSR) ####
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


#### Run sample2 ####
SAMPLENAME=sample2

#NOTE: Before go on to VCF calling,
#      run lines 49 to 138 (from fastQC to BQSR)


#-------------------------------------#
#            VCF calling              #
#-------------------------------------#
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6
module load BCFtools/1.9


### HaplotypeCaller with GVCF mode ####
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


### GenomicsDBImport ####
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
  
### GenotypeGVCFs: convert gvcf to vcf ####
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
  
  
### Merge vcfs by chromosome ####
# Extract head of vcf
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



