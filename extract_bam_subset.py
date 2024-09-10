import argparse
from bai import Bai
from Bio import bgzf
import io
import gzip
import pysam

# Block gzip end of file marker
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

# Function to initialize the BAI index and BAM header
def initialize_bam(bam_file, bai_file):
    b_idx = Bai(bai_file)  # Initialize BAI index
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        headerobject = bam.header  # Extract BAM header
    return headerobject, b_idx

# Function to read the BAM header bytes
def read_bam_header(bam_file, header_offset):
    with open(bam_file, 'rb') as f:
        header_bytes = f.read(header_offset)  # Read header bytes
    return header_bytes

# Function to query the BAM region and extract relevant blocks
def extract_bam_region(bai_index, bam_file, ref_id, start_coords, end_coords, padd):
    start = bai_index.query(ref_id, start_coords, start_coords + 1)
    end = bai_index.query(ref_id, end_coords + padd, end_coords + 1 + padd)

    start_startb, start_startoff = bgzf.split_virtual_offset(start.voffset_beg)
    end_startb, end_startoff = bgzf.split_virtual_offset(end.voffset_beg)
    end_endb, _ = bgzf.split_virtual_offset(end.voffset_end)

    with open(bam_file, 'rb') as f:
        f.seek(start_startb)
        chunk1 = f.read(end_startb - start_startb)
        chunk2 = f.read(end_endb - end_startb)

    first_block, middle_blocks = process_first_block(chunk1, start_startoff)
    last_block = process_last_block(chunk2, end_startoff)

    return first_block, middle_blocks, last_block

# Function to process the first block of the BAM region
def process_first_block(chunk1, start_offset):
    filehndl = io.BytesIO(chunk1)
    values = [x for x in bgzf.BgzfBlocks(filehndl)]
    frst_blck = chunk1[:values[0][1]]
    frst_blck_cln = gzip.compress(gzip.decompress(frst_blck)[start_offset:])
    blks_nofirst_nolast = chunk1[values[0][1]:]
    return frst_blck_cln, blks_nofirst_nolast

# Function to process the last block of the BAM region
def process_last_block(chunk2, end_offset):
    filehndl = io.BytesIO(chunk2)
    values = [x for x in bgzf.BgzfBlocks(filehndl)]
    frst_blck_end = chunk2[:values[0][1]]
    frst_blck_end_cln = gzip.compress(gzip.decompress(frst_blck_end)[:end_offset])
    return frst_blck_end_cln

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
    headerobject, b_idx = initialize_bam(bam_file, bai_file)
    
    ref_id = headerobject.get_tid(chromosome)
    print(f"Reference ID for {chromosome} is {ref_id}")

    header = b_idx.query(0, 0, 1)
    header_bytes = read_bam_header(bam_file, bgzf.split_virtual_offset(header.voffset_beg)[0])

    first_block, middle_blocks, last_block = extract_bam_region(b_idx, bam_file, ref_id, start_coords, end_coords, padd=100000)

    write_bam_region(output_file, header_bytes, first_block, middle_blocks, last_block)

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
