from BGZblocksfrombam.source.settings.log_settings import logger
import pysam


# Function to parse the region string
def parse_region(region):
    chromosome, coords = region.split(':')
    start, end = map(int, coords.split('-'))
    return chromosome, start, end


def initialize_bam(bam_file, bai_file):
    logger.info(f'input {bam_file, bai_file}')
    
    with pysam.AlignmentFile(bam_file, index_filename=bai_file) as bam:
        headerobject = bam.header  # Extract BAM header
    return headerobject


def collapse_bed(bed,padding):
    positions = []

    with open(bed) as bedfile:
        for line in bedfile:
            positions.append(line.rstrip().split('\t'))
    positions = [[x[0],int(x[1]),int(x[2])] for x in positions if len(x)==3]
    positions_request = []
    contigs = list(set([x[0] for x in positions]))

    contigs_sex = [x for x in contigs if 'y' in x.lower() or 'x' in x.lower()]
    contigs = [x for x in contigs if 'y' not in x.lower() and  'x' not in x.lower()]
    
    contigs = sorted(contigs, key=lambda s: int(s.lstrip('chr')))

    contigs = contigs+contigs_sex
    
    for contig in contigs:
        done_positions = []
        pos_on_contig = [x for x in positions if x[0]==contig]
        for n,pos in enumerate(pos_on_contig):
            if n in done_positions:
                continue
        
            # start & end pos of ROI in json
            start = pos[1]
            end = pos[2]
        
            # ROI end with padding and pot read len
            req_end = end+padding
            req_start = start-padding
            if req_start<1:
                req_start=1
        
        # searching sorted positions on contig for region included in current region 
            for nn in range(n+1,len(pos_on_contig)):
                start_next = pos_on_contig[nn][1]
                end_next = pos_on_contig[nn][2]
            
                # if start_next smaller than req_end --> new req_end is the end_next + padd ... --> regions are collapsed 
                # --> less data is requested
                start_next_req=start_next-padding
                if start_next_req<1:
                    start_next_req=1
                if start_next_req<req_end:
                    done_positions.append(nn)
                    if (end_next+padding)>req_end: 
                        req_end = end_next+padding
        
        
            positions_request.append((contig,req_start,req_end))
                
    chroms=[x[0]for x in positions_request]
    starts=[x[1]for x in positions_request]
    ends=[x[2]for x in positions_request]


    return [[x,xx,xxx] for x,xx,xxx in zip(chroms, starts, ends)]