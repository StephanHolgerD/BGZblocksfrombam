from Bio import bgzf



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