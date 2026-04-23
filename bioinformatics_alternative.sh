#dependencies
curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-*.sh
wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh
bash ./Anaconda3-2022.10-Linux-x86_64.sh
source ~/.bashrc
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install samtools
conda install bwa
conda install freebayes
conda install picard
conda install bedtools
conda install trimmomatic
conda install fastqc
sudo apt install libvcflib-tools
sudo apt install snpeff
#set up
pwd
mkdir ngs_course
mkdir ngs_course/dnaseq
cd ngs_course/dnaseq
mkdir data meta results logs
cd ~/ngs_course/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq
	#input files
wget -c https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz -O NGS0001.R1.fastq.gz
wget -c https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz -O NGS0001.R2.fastq.gz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
#set up
mv *fastq.gz ~/ngs_course/dnaseq/data/untrimmed_fastq
mv hg19.fa.gz ~/ngs_course/dnaseq/data

#pre-alignment QC
#Quality asessment and trimming
cd ~/ngs_course/dnaseq/data/untrimmed_fastq
fastqc -t 4 *.fastq.gz
mkdir ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads
mv *fastqc* ~/ngs_course/dnaseq/results/fastqc_untrimmed_reads/
trimmomatic PE \
  -threads 4 \
  -phred33 \
  ~/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
  ~/ngs_course/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
  -baseout ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.40-hdfd78af_0/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

#basic quality assessment of paired trimmed sequencing data
cd ~/ngs_course/dnaseq/data/trimmed_fastq
fastqc *P

#alignment
#Align the paired trimmed fastq files using bwa mem and reference genome hg19 
mkdir -p ~/ngs_course/dnaseq/data/reference
mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/
bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
mkdir ~/ngs_course/dnaseq/data/aligned_data
bwa mem -t 4 -v 1 \
  -R '@RG\tID:NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz \
  ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_1P \
  ~/ngs_course/dnaseq/data/trimmed_fastq/NGS0001_trimmed_2P \
  > ~/ngs_course/dnaseq/data/aligned_data/NGS0001.sam
cd ~/ngs_course/dnaseq/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam
samtools sort NGS0001.bam > NGS0001_sorted.bam
samtools index NGS0001_sorted.bam
ls

#Perform duplicate marking
picard MarkDuplicates I= NGS0001_sorted.bam O= NGS0001_Sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_Sorted_marked.bam

#Quality Filter the duplicate marked BAM file 
samtools view -F 1796  -q 20 -o NGS0001_Sorted_filtered.bam NGS0001_Sorted_marked.bam
samtools index NGS0001_Sorted_filtered.bam

#Generate standard alignment statistics 
samtools flagstat NGS0001_Sorted_filtered.bam > NGS0001_flagstats.txt
samtools view NGS0001_Sorted_filtered.bam | head -n 20
samtools idxstats NGS0001_Sorted_filtered.bam > NGS0001_idxstats.txt
picard CollectInsertSizeMetrics       I=NGS0001_Sorted_filtered.bam       O=insert_size_metrics.txt       H=insert_size_histogram.pdf       M=0.5
bedtools genomecov -ibam NGS0001_Sorted_filtered.bam > NGS0001_coverage.txt
bedtools genomecov -ibam NGS0001_Sorted_filtered.bam -g ~/ngs_course/dnaseq/data/reference/hg19.fa.gz

#Variant Calling
#freebayes
sudo apt install bcftools
bcftools mpileup -Ou -f ~/ngs_course/dnaseq/data/reference/hg19.fa   ~/ngs_course/dnaseq/data/aligned_data/NGS0001_Sorted_filtered.bam |   bcftools call -mv -Oz -o ~/ngs_course/dnaseq/results/NGS0001.vcf.gz
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001.vcf.gz

#quality Filter Variants
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
~/ngs_course/dnaseq/results/NGS0001.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf
bedtools intersect -header -wa \
  -a ~/ngs_course/dnaseq/results/NGS0001_filtered.vcf \
  -b /home/ubuntu/ngs_course/dnaseq/data/annotation.bed \
  > ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf
bgzip ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf
tabix -p vcf ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf.gz

#Variant Annotation and Prioritization
#Annotate using ANNOVAR and snpEFF
tar -zxvf annovar.latest.tar.gz
cd annovar
chmod +x annotate_variation.pl
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf.gz > ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.avinput
./table_annovar.pl ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.avinput humandb/ -buildver hg19\ 
   -out ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation -remove\ 
      -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

snpEff databases | grep -i "Human"
snpEff -Xmx4g GRCh37.75 \
   ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.vcf.gz \
   > ~/ngs_course/dnaseq/results/NGS0001_snpeff_annotated.vcf

#Perform basic variant prioritization: filter to exonic variants not seen in dbSNP 
awk -F'\t' '$6 == "exonic" && $11 == "."' ~/ngs_course/dnaseq/results/NGS0001_filtered_annotation.hg19_multianno.txt > prioritized_variants.txt
wc -l prioritized_variants.csv
