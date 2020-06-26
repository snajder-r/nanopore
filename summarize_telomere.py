import pysam
import os
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib
matplotlib.use("Agg")

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


file_dict = { 
            'germline': '../telomeres/edgecase_only_primary/Germline.tails.chopped.bam',
            'primary': '../telomeres/edgecase_only_primary/Primary.tails.chopped.bam',
            'relapse': '../telomeres/edgecase_only_primary/Relapse.tails.chopped.bam'
            }

samples = list(file_dict.keys())

all_telomeres = pd.DataFrame({
                   'read': pd.Series([], dtype='str'),
                   'sample': pd.Series([], dtype='str'),
                   'contig': pd.Series([], dtype='str'),
                   'pq': pd.Series([], dtype='str'),
                   'length': pd.Series([], dtype=int)})

for key in samples:
    with pysam.AlignmentFile(file_dict[key], "rb") as f:
        for line in f.fetch():
            is_q = line.flag & 0x8000
            all_telomeres = all_telomeres.append(
                    {
                        'read': line.query_name,
                        'sample': key,
                        'contig': line.reference_name,
                        'pq': 'q' if is_q else 'p',
                        'length': len(line.query_sequence)
                    }
                    , ignore_index=True
                )

means = all_telomeres.groupby(['read', 'sample', 'contig', 'pq']).mean()
all_telomeres.to_csv('../telomeres/edgecase/chopped_telomeres.tsv', sep='\t')
means.to_csv('../telomeres/edgecase/chopped_telomeres_mean.tsv', sep='\t')

#all_telomeres = all_telomeres.groupby(['read', 'sample','contig','pq']).sum().reset_index()
print(all_telomeres)
with PdfPages('../telomeres/edgecase/telomeres_lengths.pdf') as pdf:
    for contig in set(all_telomeres['contig']):
        for pq in ['p','q']:
            plt.figure(figsize=(7,7))
            idx = (all_telomeres['contig']==contig) & (all_telomeres['pq'] == pq)
            for i,sample in enumerate(samples):
                s_lengths = all_telomeres.loc[idx & (all_telomeres['sample']==sample)]['length']
                plt.scatter(i*2+np.random.rand(len(s_lengths)), s_lengths)
            plt.xticks(np.arange(len(samples))*2+0.5, samples, rotation=20)
            test = (scipy.stats.ttest_ind(all_telomeres.loc[idx & (all_telomeres['sample']=='germline')]['length'], all_telomeres.loc[idx & (all_telomeres['sample']=='primary')]['length']))
            plt.title('%s %s ' % (contig,pq))
            pdf.savefig()
            plt.close()

plt.figure(figsize=(8,8))
sample_telomere_dist = {sample:np.log10(all_telomeres.loc[(all_telomeres['sample']==sample)]['length']) for sample in samples}
for i,sample in enumerate(samples):
    s_lengths = sample_telomere_dist[sample]
    plt.scatter(i*2+np.random.rand(len(s_lengths)), s_lengths)
plt.xticks(np.arange(len(samples))*2+0.5, samples, rotation=20)
plt.title('Overall overhang length')
plt.savefig('../telomeres/edgecase/telomeres_lengths_all.pdf')
plt.close()

plt.title('Overall overhang length distributions (all chromosomes)')
plt.boxplot([sample_telomere_dist[sample] for sample in samples])
plt.xticks(np.arange(1,1+len(samples)), samples, rotation=20)
plt.ylabel('Length (log10)')
plt.savefig('../telomeres/edgecase/telomeres_lengths_boxplot.pdf')
print(scipy.stats.ttest_ind(sample_telomere_dist['germline'], sample_telomere_dist['primary']))
print(scipy.stats.ttest_ind(sample_telomere_dist['germline'], sample_telomere_dist['relapse']))
print(scipy.stats.ttest_ind(sample_telomere_dist['primary'], sample_telomere_dist['relapse']))

