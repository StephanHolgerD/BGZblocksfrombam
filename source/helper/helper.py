from settings.log_settings import logger
import pysam


# Function to parse the region string
def parse_region(region):
    chromosome, coords = region.split(':')
    start, end = map(int, coords.split('-'))
    return chromosome, start, end


def initialize_bam(bam_file, bai_file):
    logger.info(f'input {bam_file, bai_file}')
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        headerobject = bam.header  # Extract BAM header
    return headerobject


