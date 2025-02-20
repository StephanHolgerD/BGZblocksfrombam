
import zlib
import struct
from BGZblocksfrombam.source.bgzf.bgzf_marker import _bgzf_header

from BGZblocksfrombam.source.settings.log_settings import logger

def bgzip_block(block):
    logger.info(f'clean data block')
    
    c = zlib.compressobj(
    6, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
    )
    compressed = c.compress(block) + c.flush()
    bsize = struct.pack("<H", len(compressed) + 25)  # includes -1
    crc = struct.pack("<I", zlib.crc32(block) & 0xFFFFFFFF)
    uncompressed_length = struct.pack("<I", len(block))
    
    
    
    clean_block =  _bgzf_header + bsize + compressed + crc + uncompressed_length
    return clean_block