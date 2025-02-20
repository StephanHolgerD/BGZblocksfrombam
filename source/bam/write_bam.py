# Function to write the extracted region to a new BAM file
from BGZblocksfrombam.source.bgzf.bgzf_marker import _bgzf_eof
from BGZblocksfrombam.source.settings.log_settings import logger
import pysam
import os

def write_bam_region(output_file, header_bytes, blocks):
    
    
    with open(output_file, 'wb') as o:
        o.write(header_bytes)
        for block in blocks:
            first_block, middle_blocks, last_block = block
            o.write(first_block)
            o.write(middle_blocks)
            o.write(last_block)
        o.write(_bgzf_eof)
    logger.info(f"BAM region written to {output_file}")
    pysam.index(output_file)
    if os.path.exists(f'{output_file}.bai'):
        logger.info(f"Index written to {output_file}.bai")
    
