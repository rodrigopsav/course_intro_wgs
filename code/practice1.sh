#-------------------------------------#
#     Practice 1: file formats        #
#-------------------------------------#


### Searching for HPCC modules ####
module spider samtools
module spider SAMtools/1.10


### Loading HPCC modules ####
module purge
module load GCC/8.3.0
module load HTSlib/1.10.2
module load SAMtools/1.10
module load BCFtools/1.10


### FASTA format ####
cd
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

# grep options:
# -m 1, --max-count, stop reading a file after 1 matching line   


### Indexing FASTA with samtools ####
samtools faidx genome.fa
cat genome.fa.fai

  
### FASTQ format ####
zcat example9_S2_R1_001.R1.fastq.gz | less

  
### SAM / BAM formats ####
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

  
### VCF format ####
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

  
### GTF format ####
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

# Searching with grep
grep -v "^#" Bos_taurus.ARS-UCD1.2.104.gtf | grep -i "dgat" | awk '$3 == "gene" {print}'






