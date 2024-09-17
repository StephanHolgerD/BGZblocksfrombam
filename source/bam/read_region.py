import requests
from requests.adapters import HTTPAdapter, Retry
import gzip
from Bio import bgzf
from bam.calculate_region import get_end,get_start
from bgzf.block import bgzip_block
from bai.baiparser import get_bai_bins
import io


def read_bam_region(bai_file, bam_file, ref_id, start_coords, end_coords, padd):
    bai_bins = get_bai_bins(bai_file,ref_id)
    print(bai_bins)
    start_startb,start_startoff = get_start(bai_bins,start_coords)
    end_startb,end_startoff = get_end(bai_bins,end_coords,start_startb)
    
    
    
    if bam_file.startswith('https'):
        headers1 = {"Range": f"bytes={start_startb}-{end_startb-1}"}
        headers2 = {"Range": f"bytes={end_startb}-{end_startb+20_000}"}
        with requests.Session() as s:
            s.mount('https://', HTTPAdapter(max_retries=5))
            response1 = s.get(bam_file,headers=headers1)
            chunk1=response1.content
            response2 = s.get(bam_file,headers=headers2)
            chunk2=response2.content
    else:
        with open(bam_file,'rb')as f:
            f.seek(start_startb)
            chunk1 = f.read(end_startb-start_startb)
            chunk2 = f.read(20_000)

    filehndl = io.BytesIO(chunk1)
    values = [x for x  in bgzf.BgzfBlocks(filehndl)]
    frst_blck = chunk1[:values[0][1]]
    block = gzip.decompress(frst_blck)[start_startoff:]

    frst_blck_cln = bgzip_block(block)
    
    blks_nofirst_nolast = chunk1[values[0][1]:]



    filehndl = io.BytesIO(chunk2)

    values = []


    for n,x in enumerate(bgzf.BgzfBlocks(filehndl)):
        values.append(x)
        break

    frst_blck_end = chunk2[:values[0][1]]
    block = gzip.decompress(frst_blck_end)[:end_startoff]
    frst_blck_end_cln = bgzip_block(block)
    
    return frst_blck_cln,blks_nofirst_nolast,frst_blck_end_cln
    #blks_nofirst_nolast = chunk1[values[0][1]:]

