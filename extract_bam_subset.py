import argparse
# from bai import Bai
from Bio import bgzf
import io
import gzip
import pysam
import requests
from requests.adapters import HTTPAdapter, Retry

import zlib
import struct
from bai.baiparser import get_bai_bins, get_header_bytes

# Block gzip end of file marker
_bgzf_magic = b"\x1f\x8b\x08\x04"
_bgzf_header = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_bytes_BC = b"BC"

# Function to initialize the BAI index and BAM header
def initialize_bam(bam_file, bai_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        headerobject = bam.header  # Extract BAM header
    return headerobject

# Function to read the BAM header bytes
def read_bam_header(bam_file, header_bytes):
    if bam_file.startswith('https'):
        
        if header_bytes[0]==0:
            headers = {"Range": f"bytes=0-{20_000}"}
            print(headers)
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
            
            
            block= gzip.decompress(header_block_cln)[:header_bytes[1]]
            
            header = bgzip_block(block)
            
            
            
            
        else:
            header_bytes=header_bytes[0]
            headers = {"Range": f"bytes=0-{header_bytes-1}"}
            print(headers)
            with requests.Session() as s:
                s.mount('https://', HTTPAdapter(max_retries=5))
                response = s.get(bam_file,headers=headers)
                header = response.content
        
    else:
        with open(bam_file, 'rb') as f:
            header = f.read(header_bytes[0])  # Read header bytes
    return header

# Function to query the BAM region and extract relevant blocks


def bgzip_block(block):
    
    c = zlib.compressobj(
    6, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
    )
    compressed = c.compress(block) + c.flush()
    bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
    crc = struct.pack("<I", zlib.crc32(block) & 0xFFFFFFFF)
    uncompressed_length = struct.pack("<I", len(block))
    
    
    
    clean_block =  _bgzf_header + bsize + compressed + crc + uncompressed_length
    return clean_block
def extract_bam_region(bai_file, bam_file, ref_id, start_coords, end_coords, padd):
    bai_bins = get_bai_bins(bai_file,ref_id)
    print(bai_bins)
    start_startb,start_startoff = get_start(bai_bins,start_coords)
    end_startb,end_startoff = get_end(bai_bins,end_coords,start_startb)
    
    
    
    if bam_file.startswith('https'):
        headers1 = {"Range": f"bytes={start_startb}-{end_startb-1}"}
        headers2 = {"Range": f"bytes={end_startb}-{end_startb+20_000}"}
        
        print(headers1)
        print(headers2)
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

def get_start(bai_bins,start_coords):
    start_offset={k:v for n,(k,v) in enumerate(zip(bai_bins.keys(),bai_bins.values())) 
                  if k<=start_coords and list(bai_bins.keys())[n+1]>start_coords}
    start_startb,start_startoff = bgzf.split_virtual_offset(list(start_offset.values())[0])
    return start_startb,start_startoff

def get_end(bai_bins,end_coords,start_startb):
    end_offset_k=[k for n,(k,v) in enumerate(zip(bai_bins.keys(),bai_bins.values())) 
                  if k>end_coords and bgzf.split_virtual_offset(v)[0]>start_startb][0]
    end_offset={end_offset_k:bai_bins[end_offset_k]}
    end_startb,end_startoff = bgzf.split_virtual_offset(list(end_offset.values())[0])
    return end_startb,end_startoff

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

    first_block, middle_blocks, last_block = extract_bam_region(bai_file, bam_file, ref_id, start_coords, end_coords, padd=100000)

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
