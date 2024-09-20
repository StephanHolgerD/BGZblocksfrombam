# Function to write the extracted region to a new BAM file
from bgzf.bgzf_marker import _bgzf_eof
from settings.log_settings import logger



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
