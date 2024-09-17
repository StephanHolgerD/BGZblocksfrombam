import argparse
# from bai import Bai
from Bio import bgzf
import io
import gzip
import pysam
import requests
from requests.adapters import HTTPAdapter, Retry

from bai.baiparser import get_bai_bins, get_header_bytes
from bam.read_header import read_bam_header
from bam.read_region import read_bam_region
from bgzf.bgzf_marker import _bgzf_eof


# Function to initialize the BAI index and BAM header
def initialize_bam(bam_file, bai_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        headerobject = bam.header  # Extract BAM header
    return headerobject

# Function to write the extracted region to a new BAM file
def write_bam_region(output_file, header_bytes, first_block, middle_blocks, last_block):
    with open(output_file, 'wb') as o:
        o.write(header_bytes)
        o.write(first_block)
        o.write(middle_blocks)
        o.write(last_block)
        o.write(_bgzf_eof)
    print(f"BAM region written to {output_file}")




# Main function to handle the extraction process
def extract_bam_region_to_file(bam_file, bai_file, chromosome, start_coords, end_coords, output_file):
    headerobject = initialize_bam(bam_file, bai_file)
    
    ref_id = headerobject.get_tid(chromosome)
    print(f"Reference ID for {chromosome} is {ref_id}")

    header_end = get_header_bytes(bai_file)
    
    
    header = read_bam_header(bam_file, header_end)

    first_block, middle_blocks, last_block = read_bam_region(bai_file, bam_file, ref_id, start_coords, end_coords, padd=100000)

    write_bam_region(output_file, header, first_block, middle_blocks, last_block)

# Parse command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description='Extract a specific region from a BAM file.')
    parser.add_argument('bam_file', type=str, help='Path to the BAM file')
    parser.add_argument('bai_file', type=str, help='Path to the BAI index file')
    parser.add_argument('region', type=str, help='Region in the format "chromosome:start-end"')
    parser.add_argument('output_file', type=str, help='Path to the output BAM file')
    return parser.parse_args()

# Function to parse the region string
def parse_region(region):
    chromosome, coords = region.split(':')
    start, end = map(int, coords.split('-'))
    return chromosome, start, end

# Main entry point
if __name__ == "__main__":
    args = parse_args()
    bam_file = args.bam_file
    bai_file = args.bai_file
    region = args.region
    output_file = args.output_file

    # Parse region into chromosome, start, and end coordinates
    chromosome, start_coords, end_coords = parse_region(region)

    # Run the extraction process
    extract_bam_region_to_file(bam_file, bai_file, chromosome, start_coords, end_coords, output_file)
