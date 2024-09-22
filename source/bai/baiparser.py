from io import BytesIO
import os

from settings.log_settings import logger


def _read_int(x: BytesIO, n_bytes: int, signed: bool=True) -> int:
    return int.from_bytes(x.read(n_bytes), byteorder="little", signed=signed)

def read_uint8(x: BytesIO) -> int:
    return _read_int(x, 1, False)

def read_int8(x: BytesIO) -> int:
    return _read_int(x, 1, True)

def read_uint16(x: BytesIO) -> int:
    return _read_int(x, 2, False)

def read_int16(x: BytesIO) -> int:
    return _read_int(x, 2, True)

def read_uint32(x: BytesIO) -> int:
    return _read_int(x, 4, False)

def read_int32(x: BytesIO) -> int:
    return _read_int(x, 4, True)

def read_int64(x: BytesIO) -> int:
    return _read_int(x, 8, True)

def read_uint64(x: BytesIO) -> int:
    return _read_int(x, 8, False)

def make_virtual_offset(coffset, uoffset):
    return (coffset << 16) | uoffset
    pass

def split_virtual_offset(virtual_offset):
    coffset = virtual_offset >> 16
    uoffset = virtual_offset ^ (coffset << 16)
    return coffset, uoffset

def get_bai_bins(filename,refid):
    linear_bin_dict_refid = {}
    lin_index_binlen = 16384
    total_size = os.stat(filename).st_size
    with open(filename, "rb") as bai:
        magic = read_int32(bai)
        n_ref = read_uint32(bai)
        for i in range(n_ref):

            n_bins = read_uint32(bai) 
            for j in range(n_bins):
                bin = read_uint32(bai)
                n_chunks = read_uint32(bai)
                for k in range(n_chunks):
                    chunk_beg = read_uint64(bai)
                    chunk_end = read_uint64(bai)
                    assert make_virtual_offset(*split_virtual_offset(chunk_beg)) == chunk_beg
                    assert make_virtual_offset(*split_virtual_offset(chunk_end)) == chunk_end
            n_intv = read_uint32(bai)
            for k in range(n_intv):
                ioffset = read_uint64(bai)
                if i==refid:
                    if ioffset!=0:
                        linear_bin_dict_refid[k*lin_index_binlen]=ioffset
    return linear_bin_dict_refid   

             
def get_header_bytes(filename):
    logger.info(f'input {filename}')
    
    refid = 0
    with open(filename, "rb") as bai:
        magic = read_int32(bai)
        n_ref = read_uint32(bai)
        for i in range(n_ref):

            n_bins = read_uint32(bai) 
            for j in range(n_bins):
                bin = read_uint32(bai)
                n_chunks = read_uint32(bai)
                for k in range(n_chunks):
                    chunk_beg = read_uint64(bai)
                    chunk_end = read_uint64(bai)
                    assert make_virtual_offset(*split_virtual_offset(chunk_beg)) == chunk_beg
                    assert make_virtual_offset(*split_virtual_offset(chunk_end)) == chunk_end
            n_intv = read_uint32(bai)
            
            for k in range(n_intv):
                ioffset = read_uint64(bai)
                return split_virtual_offset(ioffset)