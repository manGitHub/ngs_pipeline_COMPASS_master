# usage:
# python3 cnvkit_QCI.py \
#/data/Compass/dev/valdezkm/ngs/CP12446/MD-23-2205/CP12446_T1D_E/cnvkit/CP12446_T1D_E.call.cns \
#/data/Compass/Ref/hg19/design/NextraExome_v1.target.hg19.merged.bed \
#/data/Compass/dev/NextraExome/CNVPanel/FinalExomeCNVgenelist_459_081224.xlsx \
#/data/Compass/dev/valdezkm/cnvkit/QCI.blank.CN.vcf \
#CP12446_T1D_E \
#CP12446_T1D_E.cnvkit.vcf

import pandas as pd
import sys
import numpy as np
import os

cns_file = sys.argv[1]
bed_file = sys.argv[2]
gene_file = sys.argv[3]
header_file = sys.argv[4]
samp = sys.argv[5]
out_file = sys.argv[6]

if os.stat(cns_file).st_size == 0:
    open(out_file, 'w').close()
    exit()

choose_cn = lambda x: np.median(x[~np.isnan(x)])
gene_aliases = {
     "AFDN":"MLLT4",
    "NSD2":"WHSC1",
    "NSD3":"WHSC1L1",
    "SEPTIN9":"SEPT9",
    "TENT5C":"FAM46C",
    "H3-3A":"H3F3A",
    "H3C2":"HIST1H3B",
            }
reverse_alias = {v: k for k, v in gene_aliases.items()}
gene_aliases['TERC'] = np.nan

# genes
gene_df = pd.read_excel(gene_file, header = None, names = ['gene'])
gene_df = gene_df.replace(gene_aliases)
gene_df = gene_df.dropna()
genes = gene_df['gene'].to_list()

# bed file
bed = pd.read_csv(bed_file, sep='\t',header = None, names = ['chrom','start','end','gene','info'])
bed = bed.drop('gene', axis=1).join(bed['gene'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('gene'))
bed['gene'] = bed['gene'].str.replace('___.*','',regex=True)
bed = bed.sort_values(by=['chrom','start','end']).groupby('gene',sort= False,as_index=False).agg({'chrom': 'first', 'start': 'first','end': 'last'})
bed = bed[bed['gene'].isin(genes)]
bed = bed.replace(reverse_alias)

# cnvkit cns file
cns_raw = pd.read_csv(cns_file, sep='\t')
cns = cns_raw.drop('gene', axis=1).join(cns_raw['gene'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('gene'))
cns = cns[['cn','cn1','cn2','gene']]
cns['gene'] = cns['gene'].str.replace('___.*','',regex=True)
cns = cns.groupby('gene',sort= False,as_index=False).agg({'cn': 'median','cn1':choose_cn,'cn2':choose_cn})
# if cn1 + cn2 != cn, change to NA
cns.loc[cns['cn1']+cns['cn2']!=cns['cn'],['cn1','cn2']] = np.nan
cns = cns[cns['gene'].isin(genes)]
cns = cns.replace(reverse_alias)

# merge data
dat = bed.merge(cns,on='gene',how='inner')
# add CDKN2B, if segment not found then take CDKN2A
n2b_coord = bed[bed['gene'].isin(['CDKN2B'])]
n2b_cn = cns_raw[(cns_raw['chromosome'].isin([n2b_coord['chrom'].to_string(index=False).strip()])) & (cns_raw['start']<=int(n2b_coord['start'].to_string(index=False).strip())) & (cns_raw['end']>=int(n2b_coord['end'].to_string(index=False).strip()))]
if n2b_cn.empty:
    n2b = cns[cns['gene'].isin(['CDKN2A'])]
    n2b['gene'] = 'CDKN2B'
else:
    n2b_cn = n2b_cn[['cn','cn1','cn2']]
    n2b_cn = n2b_cn.groupby(level=0).median()
    n2b = pd.DataFrame(['CDKN2B',n2b_coord['chrom'].to_string(index=False).strip(),n2b_coord['start'].to_string(index=False).strip(),n2b_coord['end'].to_string(index=False).strip(),n2b_cn['cn'].to_string(index=False).strip(),n2b_cn['cn1'].to_string(index=False).strip(),n2b_cn['cn2'].to_string(index=False).strip()]).T
    n2b.columns = dat.columns
dat = pd.concat([dat,n2b])
dat['cn'] = np.floor(pd.to_numeric(dat['cn'], errors='coerce')).astype('Int64')
dat['cn1'] = np.floor(pd.to_numeric(dat['cn1'], errors='coerce')).astype('Int64')
dat['cn2'] = np.floor(pd.to_numeric(dat['cn2'], errors='coerce')).astype('Int64')

# header 
h = open(header_file,'r')
header = h.read()
out = open(out_file,'w')
out.write(header)
out.close()

# create vcf
dat['ID'] = '.'
dat['REF'] = 'N'
dat['ALT'] = '<CNV>'
dat['QUAL'] = '.'
dat['FILTER'] = 'PASS'
dat['INFO'] = 'SVTYPE=CNV;END=' + dat['end'].astype(str) + ';ANT=' + dat['gene'] + ';CN1=' + dat['cn1'].astype(str) + ';CN2=' + dat['cn2'].astype(str)
dat['FORMAT'] = 'GT:CN'
dat = dat[['chrom','start','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','cn']]
dat['cn'] = './.:' + dat['cn'].astype(str)
dat.columns = ['#CHROM','POS'] + list(dat.columns[2:-1]) + [samp]
dat.to_csv(out_file,mode='a',index=False,sep='\t')

