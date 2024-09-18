
import zlib
import struct
from bgzf.bgzf_marker import _bgzf_header

import logging
from settings import log_settings  
logging.basicConfig(
        stream=log_settings.stream,
        level=log_settings.level,
        format=log_settings.format
        )
def bgzip_block(block):
    logging.info(f'clean data block')
    
    c = zlib.compressobj(
    6, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
    )
    compressed = c.compress(block) + c.flush()
    bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
    crc = struct.pack("<I", zlib.crc32(block) & 0xFFFFFFFF)
    uncompressed_length = struct.pack("<I", len(block))
    
    
    
    clean_block =  _bgzf_header + bsize + compressed + crc + uncompressed_length
    return clean_block