# Simulate reads using wgsim
wgsim -N 500000 GCA_000001405.15_GRCh38_no_alt_analysis_set.fna testdata.read1.fq testdata.read2.fq

# Gzip the reads
gzip testdata.*

# Align using bwa and sort with samtools
bwa mem GCA_000001405.15_GRCh38_no_alt_analysis_set.fna testdata.read1.fq.gz testdata.read2.fq.gz | samtools sort -o testdata.bam

# Index the bam file
samtools index testdata.bam

# Remove the fastq files
rm testdata.read1.fq.gz testdata.read2.fq.gz