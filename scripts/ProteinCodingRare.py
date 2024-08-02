# To replace GeneAnnotation.pl
# Written by Rob Schultz 9/23
#  Takes in the temporary output from GeneAnnotation.py, Creates coding.rare annotation

import sys
import pandas as pd
import csv

blacklist_path =  sys.argv[1]
whitelist_path =  sys.argv[2]
anno_path =  sys.argv[3]
maf =  float(sys.argv[4])
tmp_anno_path =  sys.argv[5]

blacklist_header = ['Chr','Start','End','Ref','Alt','Info']
whitelist_header = ['Chr','Start','End','Gene']
black_df = pd.read_csv(blacklist_path,sep = '\t', names = blacklist_header, dtype = str) 
white_df = pd.read_csv(whitelist_path,sep = '\t', names = whitelist_header, dtype = str)
anno_df = pd.read_csv(anno_path,sep = '\t', dtype = str, keep_default_na = False,quoting = csv.QUOTE_NONE)

# pop_freqs = [i for i in anno_df.columns][12:29]
pop_freqs =['1000g2014oct_all', '1000g2014oct_eur', '1000g2014oct_afr', '1000g2014oct_amr', '1000g2014oct_eas', '1000g2014oct_sas', 'esp6500_all', 'esp6500_ea', 'esp6500_aa', 'ExAC_ALL_nonTCGA', 'ExAC_AFR_nonTCGA', 'ExAC_AMR_nonTCGA', 'ExAC_EAS_nonTCGA', 'ExAC_FIN_nonTCGA', 'ExAC_NFE_nonTCGA', 'ExAC_OTH_nonTCGA', 'ExAC_SAS_nonTCGA']

keepers = []
#anno_df[pop_freqs] = anno_df[pop_freqs].applymap(lambda s: '-2' if s=='.' else s)
#anno_df[pop_freqs] = anno_df[pop_freqs].applymap(lambda s: -2 if s=='.' else float(s))
for i in anno_df.index:
    isblacklist = False
    iswhitelist = False
    row = anno_df.loc[i]
    for b in black_df.index:
        if all(row[blacklist_header[:-1]] == black_df.loc[b,blacklist_header[:-1]]):
            isblacklist = True
            break
    if isblacklist:
        continue
    for w in white_df.index:
        #if row['Chr'] == white_df.loc[w,'Chr'] and row['Start'] >= white_df.loc[w,'Start'] and row['Start'] <= white_df.loc[w,'End']:
        if row['Chr'] == white_df.loc[w,'Chr'] and int(row['Start']) >= int(white_df.loc[w,'Start']) and int(row['Start']) <= int(white_df.loc[w,'End']):
            keepers.append(i)
            iswhitelist = True
            break
    if iswhitelist:
        continue
    if row['Func_refGene'].startswith('splicing') or (row['Func_refGene'].startswith('exonic') and any([func in row['ExonicFunc_refGene'] for func in ['nonsynonymous','startloss','startgain','stop','frameshift']])):
        #row[pop_freqs] = row[pop_freqs].apply(lambda s: -2 if s=='.' else float(s))
        #if all(row[pop_freqs].apply(lambda s: s <= maf or pd.isnull(s))):
        if all(row[pop_freqs].apply(lambda s: float(s) <= maf if not s == '.' else True)):
            keepers.append(i)
   



rare_df = anno_df.loc[keepers]
#print(len(keepers))
#print(tmp_anno_path)


rare_df.to_csv(tmp_anno_path,sep = '\t',index = False,quoting = csv.QUOTE_NONE)

