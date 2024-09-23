
* used  https://github.com/oyvinev/bam-inspection bai parser

# requirements
* pysam 
* biopython
* io
* gzip

# Script usage

```bash
usage: extract_bam_subset.py region [-h] -i BAM_FILE -x BAI_FILE -o OUTPUT_FILE -r REGION

options:
  -h, --help            show this help message and exit

required arguments:
  -i BAM_FILE, --bam_file BAM_FILE
                        Pos. sorted and indexed bam file
  -x BAI_FILE, --bai_file BAI_FILE
                        Pos. sorted and indexed bam file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output BAM file
  -r REGION, --region REGION
                        Region in the format "chromosome:start-end"

```



```bash
usage: extract_bam_subset.py bed [-h] -i BAM_FILE -x BAI_FILE -o OUTPUT_FILE -b BED

options:
  -h, --help            show this help message and exit

required arguments:
  -i BAM_FILE, --bam_file BAM_FILE
                        Pos. sorted and indexed bam file
  -x BAI_FILE, --bai_file BAI_FILE
                        Pos. sorted and indexed bam file
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output BAM file
  -b BED, --bed BED     BED file of regions for extraction


```