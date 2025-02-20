import requests
from requests.adapters import HTTPAdapter, Retry
import io
from Bio import bgzf
import gzip
from BGZblocksfrombam.source.bgzf.block import bgzip_block
from BGZblocksfrombam.source.settings.log_settings import logger


# Function to read the BAM header bytes
def read_bam_header(bam_file, header_bytes):
    logger.info(f'input {bam_file, header_bytes}')
    
    if bam_file.startswith('https'):
        logger.info(f'BAM file  {bam_file} is remote')
        
        
        if header_bytes[0]==0:
            logger.info(f'BAM header and first data block mixed')
            
            headers = {"Range": f"bytes=0-{20_000}"}
            logger.info(f'request header with {headers} to get first block which includes the header')

            with requests.Session() as s:
                s.mount('https://', HTTPAdapter(max_retries=5))
                response = s.get(bam_file,headers=headers)

            values = []
            
            header_block = response.content
            filehndl = io.BytesIO(header_block)

            for n,x in enumerate(bgzf.BgzfBlocks(filehndl)):
                values.append(x)
                break
            header_block_cln = header_block[:values[0][1]]
            
            
            logger.info(f'decompress first block and read until offset where header ends {header_bytes[1]}')
            block= gzip.decompress(header_block_cln)[:header_bytes[1]]
            
            header = bgzip_block(block)
            
            
            
            
        else:
            logger.info(f'BAM header and first data block not mixed')
            
            header_bytes=header_bytes[0]
            headers = {"Range": f"bytes=0-{header_bytes-1}"}
            logger.info(f'request header with {headers}')
            
            with requests.Session() as s:
                s.mount('https://', HTTPAdapter(max_retries=5))
                response = s.get(bam_file,headers=headers)
                header = response.content
        
    else:
        with open(bam_file, 'rb') as f:
            header = f.read(header_bytes[0])  # Read header bytes
    return header


