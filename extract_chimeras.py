import pysam
import sys
import os
import pandas as pd

class ChimericReadAlignment:
    def __init__(self, contig, is_primary, start, end):
        self.contig = contig
        self.is_primary = is_primary
        self.start = start
        self.end = end

def getChimeraMap(f):
    ret = dict()
    for read in f.fetch():
        if not read.query_name in ret.keys():
            ret[read.query_name] = list()
        is_primary = (read.flag & 2048 == 0)
        ret[read.query_name].append(
                ChimericReadAlignment(contig = read.reference_name, 
                                      is_primary = is_primary, 
                                      start = read.reference_start, 
                                      end = read.reference_end)
                )
    return ret

def loadAnchors(anchor_file):
    anchors = pd.read_csv(anchor_file, sep='\t', header=0, skiprows=1)
    return anchors

def filterChimeraOverAnchor(chimera_map, anchors):
    filtered_map = dict()
    for (k,vl) in chimera_map.items():
        for v in vl:
            mask = v.contig == anchors['rname']
            mask = mask & (v.start <= anchors['pos'])
            mask = mask & (v.end >= anchors['pos'])
            overlaps_anchor = mask.sum() > 0
            if overlaps_anchor:
                filtered_map[k] = vl
                break
    return filtered_map

def main(argv):
    indir = argv[1]
    anchor_file = argv[2]
    anchors = loadAnchors(anchor_file)

    full_chimera_map = dict()
    fl =  [fn for fn in os.listdir(indir) if fn.endswith('.bam')]
    for i, fname in enumerate(fl):
        print('%s (%d/%d)' % (fname, i, len(fl)))
        infile = os.path.join(indir, fname)
        with pysam.AlignmentFile(infile, 'rb') as f:
            part_map = getChimeraMap(f)
            for read_name in part_map.keys():
                if not read_name in full_chimera_map.keys():
                    full_chimera_map[read_name] = list()
                full_chimera_map[read_name] = full_chimera_map[read_name]  + part_map[read_name]
        if i == 10:
            break
    print(len(full_chimera_map.keys()))

    only_chimeric = {k:full_chimera_map[k] for k in full_chimera_map.keys() if len(full_chimera_map[k])>1}
    print(len(only_chimeric.keys()))

    telomere_chimeras = filterChimeraOverAnchor(only_chimeric, anchors)
    print(len(telomere_chimeras.keys()))
    print(telomere_chimeras)

if __name__=='__main__':
    main(sys.argv)
