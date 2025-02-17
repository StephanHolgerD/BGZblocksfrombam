import argparse
import sys
def parse_args():
    parser = argparse.ArgumentParser(description='extract region from bam local or remote')
    subparsers = parser.add_subparsers(help='extract region from bam file using region string or bed')

    region_parser = subparsers.add_parser("region")
    bed_parser = subparsers.add_parser("bed")

    requiredNamed = region_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i','--bam_file', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-x','--bai_file', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-o','--output_file', help='Path to the output BAM file', required=True)
    requiredNamed.add_argument('-r','--region', type=str, help='Region in the format "chromosome:start-end"', required=True)
    
    
    requiredNamed = bed_parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-i','--bam_file', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-x','--bai_file', help='Pos. sorted and indexed bam file', required=True)
    requiredNamed.add_argument('-o','--output_file', help='Path to the output BAM file', required=True)
    requiredNamed.add_argument('-b','--bed', type=str, help='BED file of regions for extraction', required=True)
    
    
    
    

    # optArguments = bed_parser.add_argument_group('optional arguments')
    # optArguments.add_argument('--threads',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)
    # optArguments.add_argument('-p','--padding', help='number of nt around the region', default=25,type=int)


    

    # optArguments = region_parser.add_argument_group('optional arguments')
    # optArguments.add_argument('--threads',default=1, help="number of cpu's  to run in paralell, ROI <1000 will always use 1 core",type=int)
    # optArguments.add_argument('-p','--padding', help='number of nt around the region', default=25,type=int)
    print(sys.argv)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    elif len(sys.argv)==2:
        if sys.argv[-1]=='region':
            region_parser.print_help(sys.stderr)
            sys.exit(1)
        if sys.argv[-1]=='bed':
            bed_parser.print_help(sys.stderr)
            sys.exit(1)
        
    else:
        return parser.parse_args()
