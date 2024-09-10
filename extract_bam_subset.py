from bai import Bai
from Bio import bgzf
import io
import gzip
import pysam

# Block gzip end of file marker
_bgzf_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"

# Function to initialize the BAI index and BAM header
def initialize_bam(bam_file, bai_file):
    """
    Initializes BAM and BAI file handlers, and retrieves the BAM header.

    Parameters:
    bam_file (str): Path to the BAM file
    bai_file (str): Path to the BAI file

    Returns:
    headerobject: BAM header object
    b_idx: Bai index object
    """
    b_idx = Bai(bai_file)  # Initialize BAI index
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        headerobject = bam.header  # Extract BAM header
    return headerobject, b_idx

# Function to read the BAM header bytes
def read_bam_header(bam_file, header_offset):
    """
    Reads the header bytes from a BAM file.

    Parameters:
    bam_file (str): Path to the BAM file
    header_offset (int): Offset for the header block

    Returns:
    bytes: Header bytes read from the BAM file
    """
    with open(bam_file, 'rb') as f:
        header_bytes = f.read(header_offset)  # Read header bytes
    return header_bytes

# Function to query the BAM region and extract relevant blocks
def extract_bam_region(bai_index, bam_file, ref_id, start_coords, end_coords, padd):
    """
    Extracts data chunks from a BAM file based on start and end coordinates.

    Parameters:
    bai_index (Bai): BAI index object
    bam_file (str): Path to the BAM file
    ref_id (int): Reference ID for the chromosome
    start_coords (int): Starting coordinate of the region
    end_coords (int): Ending coordinate of the region
    padd (int): Padding to extend the region for read extraction

    Returns:
    tuple: Processed first block, middle block, and end block
    """
    # Query the start and end regions from the index
    start = bai_index.query(ref_id, start_coords, start_coords + 1)
    end = bai_index.query(ref_id, end_coords + padd, end_coords + 1 + padd)

    # Extract virtual offsets from queries
    start_startb, start_startoff = bgzf.split_virtual_offset(start.voffset_beg)
    end_startb, end_startoff = bgzf.split_virtual_offset(end.voffset_beg)
    end_endb, _ = bgzf.split_virtual_offset(end.voffset_end)

    # Read relevant chunks from BAM file
    with open(bam_file, 'rb') as f:
        f.seek(start_startb)  # Go to the start position
        chunk1 = f.read(end_startb - start_startb)  # Read chunk between start and end
        chunk2 = f.read(end_endb - end_startb)  # Read chunk after end

    # Process the first and last blocks
    first_block, middle_blocks = process_first_block(chunk1, start_startoff)
    last_block = process_last_block(chunk2, end_startoff)

    return first_block, middle_blocks, last_block

# Function to process the first block of the BAM region
def process_first_block(chunk1, start_offset):
    """
    Processes the first block of the extracted BAM region.

    Parameters:
    chunk1 (bytes): First chunk of BAM data
    start_offset (int): Offset within the first block

    Returns:
    tuple: Cleaned first block and remaining blocks between first and last
    """
    filehndl = io.BytesIO(chunk1)
    values = [x for x in bgzf.BgzfBlocks(filehndl)]
    frst_blck = chunk1[:values[0][1]]  # Get the first block
    frst_blck_cln = gzip.compress(gzip.decompress(frst_blck)[start_offset:])  # Clean block with offset
    blks_nofirst_nolast = chunk1[values[0][1]:]  # Remaining blocks between first and last
    return frst_blck_cln, blks_nofirst_nolast

# Function to process the last block of the BAM region
def process_last_block(chunk2, end_offset):
    """
    Processes the last block of the extracted BAM region.

    Parameters:
    chunk2 (bytes): Last chunk of BAM data
    end_offset (int): Offset within the last block

    Returns:
    bytes: Cleaned last block
    """
    filehndl = io.BytesIO(chunk2)
    values = [x for x in bgzf.BgzfBlocks(filehndl)]
    frst_blck_end = chunk2[:values[0][1]]  # Get the first block from chunk2
    frst_blck_end_cln = gzip.compress(gzip.decompress(frst_blck_end)[:end_offset])  # Clean block with offset
    return frst_blck_end_cln

# Function to write the extracted region to a new BAM file
def write_bam_region(output_file, header_bytes, first_block, middle_blocks, last_block):
    """
    Writes the extracted BAM region to a new BAM file.

    Parameters:
    output_file (str): Path to the output BAM file
    header_bytes (bytes): Header bytes of the BAM file
    first_block (bytes): First block of the BAM region
    middle_blocks (bytes): Middle blocks between the first and last
    last_block (bytes): Last block of the BAM region
    """
    with open(output_file, 'wb') as o:
        o.write(header_bytes)  # Write header
        o.write(first_block)  # Write first block
        o.write(middle_blocks)  # Write middle blocks
        o.write(last_block)  # Write last block
        o.write(_bgzf_eof)  # Write EOF marker
    print(f"BAM region written to {output_file}")

# Main function to handle the extraction process
def extract_bam_region_to_file(bam_file, bai_file, chromosome, start_coords, end_coords, output_file):
    """
    Main function to extract a specific BAM region and write it to a new BAM file.

    Parameters:
    bam_file (str): Path to the BAM file
    bai_file (str): Path to the BAI file
    chromosome (str): Chromosome name
    start_coords (int): Starting coordinate of the region
    end_coords (int): Ending coordinate of the region
    output_file (str): Path to the output BAM file
    """
    # Initialize BAM and BAI files
    headerobject, b_idx = initialize_bam(bam_file, bai_file)
    
    # Get reference ID for the chromosome
    ref_id = headerobject.get_tid(chromosome)
    print(f"Reference ID for {chromosome} is {ref_id}")

    # Read BAM header bytes
    header = b_idx.query(0, 0, 1)
    header_bytes = read_bam_header(bam_file, bgzf.split_virtual_offset(header.voffset_beg)[0])

    # Extract the BAM region
    first_block, middle_blocks, last_block = extract_bam_region(b_idx, bam_file, ref_id, start_coords, end_coords, padd=100000)

    # Write the extracted region to a new BAM file
    write_bam_region(output_file, header_bytes, first_block, middle_blocks, last_block)

# Example usage
if __name__ == "__main__":
    bam_file = '../../download/LB24-ONTCCMJH298-ready_5a8df838-5170-406f-b305-5a402ee54014.bam'
    bai_file = '../../download/LB24-ONTCCMJH298-ready_5a8df838-5170-406f-b305-5a402ee54014.bam.bai'
    chromosome = 'chr4'
    start_coords = 5000000
    end_coords = 10000000
    output_file = 'ROI.bam'

    extract_bam_region_to_file(bam_file, bai_file, chromosome, start_coords, end_coords, output_file)
