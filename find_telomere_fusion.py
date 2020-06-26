import pysam
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


file_dict = { 
            'germline': '../telomeres/edgecase/Germline.tails.merged.bam',
            'primary': '../telomeres/edgecase/Primary.tails.merged.bam',
            'relapse': '../telomeres/edgecase/Relapse.tails.merged.bam'
            }

samples = list(file_dict.keys())

fusions = pd.DataFrame({'contig_a':pd.Series([], dtype='str'),
                        'contig_b':pd.Series([], dtype='str'),
                        'count': pd.Series([], dtype=int)})

chrom_reads_good_qual = dict()
chrom_reads_bad_qual = dict()

bad_q_reads = {}
non_chimeric_reads = {}
chimeric_reads = {}
fusion_candidates = {}
fusion_candidates_unique = {}

def same_chrom(a,b):
    if a == b or a in b or b in a:
        return True
    return False

for key in samples:
    bad_q_reads[key] = 0
    non_chimeric_reads[key] = 0
    chimeric_reads[key] = 0
    fusion_candidates[key] = 0
    fusion_candidates_unique[key] = 0
    with pysam.AlignmentFile(file_dict[key], "rb") as f:
        for line in f.fetch(until_eof=True):
            if (line.flag & 0x0800)!=0: #supplementary
                continue

            if line.mapping_quality < 60:
                bad_q_reads[key] += 1
                if not line.reference_name in chrom_reads_bad_qual.keys():
                    chrom_reads_bad_qual[line.reference_name] = 0
                chrom_reads_bad_qual[line.reference_name] +=1
                continue

            if not line.reference_name in chrom_reads_good_qual.keys():
                chrom_reads_good_qual[line.reference_name] = 0
            chrom_reads_good_qual[line.reference_name] +=1

            if not line.has_tag('SA'):
                non_chimeric_reads[key] += 1
                continue
            sa = line.get_tag('SA')
            if sa == '':
                non_chimeric_reads[key] += 1
                continue
            chimeric_reads[key] += 1
            sa = sa.split(';')
            line_fusion_candidates = 0
            last_candidate_contig = None
            for sa_part in sa:
                if sa_part == '':
                    continue
                sa_part = sa_part.split(',')
                other_contig = sa_part[0]
                other_pos = sa_part[1]
                other_strand = sa_part[2]
                other_cigar = sa_part[3]
                mapQ = float(sa_part[4])
                if mapQ >= 60 and not same_chrom(other_contig, line.reference_name):
                    line_fusion_candidates += 1
                    last_candidate_contig = other_contig
            if line_fusion_candidates >0:
                fusion_candidates[key]+=1
            if line_fusion_candidates == 1:
                fusion_candidates_unique[key]+=1
                sorted_candidates = sorted([line.reference_name, last_candidate_contig])
                fusions = fusions.append({'contig_a':sorted_candidates[0], 'contig_b':sorted_candidates[1], 'count':1},ignore_index=True)


    print(key)
    print(bad_q_reads, non_chimeric_reads, chimeric_reads, fusion_candidates, fusion_candidates_unique)

    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(fusions.groupby(['contig_a', 'contig_b']).sum())


bars = []
ind = np.arange(len(samples))
#p1 = plt.bar(ind, [bad_q_reads[key] for key in samples])
p2 = plt.bar(ind, [non_chimeric_reads[key] for key in samples])
p3 = plt.bar(ind, [chimeric_reads[key] for key in samples])
p4 = plt.bar(ind, [fusion_candidates[key] for key in samples])
p5 = plt.bar(ind, [fusion_candidates_unique[key] for key in samples])
plt.legend((p2,p3,p4,p5), ('Non-chimeric alignments', 'Same-chromosome chimeric alignments', 'Multi-chromosome chimeric alignments', 'Potential fusion alignment' ))
plt.xticks(ind, samples)
plt.savefig('../telomeres/edgecase/supplementary_alignments_overview.pdf')
plt.close()

plt.title('Reads in telomere per chromosome')
chroms = sorted(list(set(chrom_reads_good_qual.keys()).union(set(chrom_reads_bad_qual.keys()))))
ind = np.arange(len(chroms))
plt.bar(ind, [chrom_reads_good_qual[c] if c in chrom_reads_good_qual.keys() else 0 for c in chroms])
#p2 = plt.bar(ind, [chrom_reads_bad_qual[c] if c in chrom_reads_bad_qual.keys() else 0 for c in chroms])
plt.xticks(ind, chroms,rotation=90)
#plt.legend((p1,p2), ('Pass QC', 'Fail QC'))
plt.tight_layout()
plt.savefig('../telomeres/edgecase/chrom_numbers.pdf')
plt.close()
