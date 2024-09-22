import pysam
import subprocess

chrom='chr2'
start=1000000
end=10000000


test_bam = 'data/testdata.bam'
test_bai = 'data/testdata.bam.bai'

pysamids = set()
with pysam.AlignmentFile('data/testdata.bam') as pysambam:
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
        
        
        
assert len(pysamids)==len(pysamids.intersection(range_request_ids))