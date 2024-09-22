import pysam
import subprocess

chrom='2'
start=1000000
end=10000000


test_bam = 'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/high_coverage_alignment/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam'
test_bai = 'data/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam.bai'

pysamids = set()
with pysam.AlignmentFile(test_bam,index_filename=test_bai) as pysambam:
    for alignment in pysambam.fetch(chrom,start,end):
        pysamids.add(alignment.query_name)
        
        




bashCommand = f'python source/extract_bam_subset.py region -i {test_bam} -x {test_bai} -r {chrom}:{start}-{end} -o test/test.bam'
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
output, error = process.communicate()
# print(output)
# print(error)


bashCommand = f'samtools index test/test.bam'
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
output, error = process.communicate()
# print(output)
# print(error)



range_request_ids = set()
with pysam.AlignmentFile('test/test.bam') as pysambam:
    for alignment in pysambam.fetch(chrom,start,end):
        range_request_ids.add(alignment.query_name)
        
        
print(len(pysamids))
print(len(range_request_ids))       
       
assert len(pysamids)==len(pysamids.intersection(range_request_ids))