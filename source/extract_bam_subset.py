import argparse
# from bai import Bai
from Bio import bgzf
import io
import gzip
import pysam
import requests
from requests.adapters import HTTPAdapter, Retry
from helper.helper import collapse_bed
from bai.baiparser import get_bai_bins, get_header_bytes
from bam.read_header import read_bam_header
from bam.read_region import read_bam_region
from bgzf.bgzf_marker import _bgzf_eof



from multiprocessing import Pool

import tempfile 

from arg_parser.arg_parser import parse_args
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


def write_bam_regions(output_file, header_bytes, blocks):
    with open(output_file, 'wb') as o:
        o.write(header_bytes)
        for x in blocks:
            first_block, middle_blocks, last_block = x
            o.write(first_block)
            o.write(middle_blocks)
            o.write(last_block)
        o.write(_bgzf_eof)
    print(f"BAM region written to {output_file}")

def get_remote_bai(bai_file):
    x = requests.get(bai_file)
    temp = tempfile.NamedTemporaryFile() 

    with open(temp.name,'wb') as out:
        out.write(x.content)
    return temp
# Main function to handle the extraction process
def extract_bam_region_to_file(bam_file, bai_file, chromosome, start_coords, end_coords, output_file):
    headerobject = initialize_bam(bam_file, bai_file)
    ref_id = headerobject.get_tid(chromosome)
    print(f"Reference ID for {chromosome} is {ref_id}")

    if 'https' in bai_file:

        temp=get_remote_bai(bai_file)
        bai_file=temp.name
        print(f"bai file  is remote, create {bai_file}")


    header_end = get_header_bytes(bai_file)
    
    
    header = read_bam_header(bam_file, header_end)

    first_block, middle_blocks, last_block = read_bam_region(bai_file, bam_file, ref_id, start_coords, end_coords)

    write_bam_region(output_file, header, first_block, middle_blocks, last_block)


def extract_bam_bed_to_file(bam_file, bai_file,bedfile, output_file,padding=5000,processes=20):
    blocks = []
    multi_args = []

    headerobject = initialize_bam(bam_file, bai_file)



    if 'https' in bai_file:

        temp=get_remote_bai(bai_file)
        bai_file=temp.name
        print(f"bai file  is remote, create {bai_file}")


    header_end = get_header_bytes(bai_file)
    header = read_bam_header(bam_file, header_end)

    regions =  collapse_bed(bedfile,padding)

    for reg in regions:

        chromosome, start_coords, end_coords = reg
        start_coords = int(start_coords)
        end_coords = int(end_coords)
        ref_id = headerobject.get_tid(chromosome)
        print(f"Reference ID for {chromosome} is {ref_id}")


        multi_args.append((bai_file, bam_file, ref_id, start_coords, end_coords))

    with Pool(processes=processes) as p:

        blocks = p.starmap(read_bam_region,multi_args)

        #first_block, middle_blocks, last_block = read_bam_region(bai_file, bam_file, ref_id, start_coords, end_coords)
        #blocks.append((first_block, middle_blocks, last_block))


    temp.close()
    temp = tempfile.NamedTemporaryFile() 

    write_bam_regions(temp.name, header, blocks)
    pysam.sort("-o", output_file, temp.name)
    pysam.index(output_file)
    temp.close()





# # Parse command-line arguments
# def parse_args():
#     parser = argparse.ArgumentParser(description='Extract a specific region from a BAM file.')
#     parser.add_argument('bam_file', type=str, help='Path to the BAM file')
#     parser.add_argument('bai_file', type=str, help='Path to the BAI index file')
#     parser.add_argument('region', type=str, help='Region in the format "chromosome:start-end"')
#     parser.add_argument('output_file', type=str, help='Path to the output BAM file')
#     return parser.parse_args()

# Function to parse the region string
def parse_region(region):
    chromosome, coords = region.split(':')
    start, end = map(int, coords.split('-'))
    return chromosome, start, end

import sys
# Main entry point
if __name__ == "__main__":
    args = parse_args()
    if args.__contains__('regions'):
        bam_file = args.bam_file
        bai_file = args.bai_file
        region = args.region
        output_file = args.output_file

    # Parse region into chromosome, start, and end coordinates
        chromosome, start_coords, end_coords = parse_region(region)

    # Run the extraction process
        extract_bam_region_to_file(bam_file, bai_file, chromosome, start_coords, end_coords, output_file)

    if args.__contains__('bed'):
        bam_file = args.bam_file
        bai_file = args.bai_file
        bed = args.bed
        output_file = args.output_file

    # Parse region into chromosome, start, and end coordinates
        # chromosome, start_coords, end_coords = parse_region(region)

    # Run the extraction process
        extract_bam_bed_to_file(bam_file, bai_file, bed, output_file)
    