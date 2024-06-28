# To replace makeDBVariantFile.pl
# Written by Rob Schultz 10/23
# creates a temporary vcf like annotation output that contains no header, and each row of the
# input with caller, qual, fs, and reads info at the end.
# 
# From perl script:
# This script takes the Patient ID and make a table of all the variants in the following format to be loaded to the database
# Chr\tStart\tEnd\tRef\tAlt\tCaller\tQUAL\tFS\tSample\tTotalReads\tAltReads
# For every library not patient
#  Caller, QUAL and FS are comma seperated
# Will take a list of file in following order.
# muTect expressed or mutect
# Strelka Expressed or Strelka (SNPS)
# Strelka Indels
# Haplotype Caller
# Platypus
# MPG
# Haplotype Caller RNASeq
# If a variant is found in first 3 files it should be rejected from next 3 files. (Somatic)
# ./makeDBVariantFile.pl NCI0002/Sample_NCI0002tumor_E_C4PJ0ANXX/calls/Sample_NCI0002tumor_E_C4PJ0ANXX.MuTect.annotated.expressed.txt NCI0002/Sample_NCI0002tumor_E_C4PJ0ANXX/calls/Sample_NCI0002tumor_E_C4PJ0ANXX.strelka.snvs.annotated.expressed.txt NCI0002/Sample_NCI0002tumor_E_C4PJ0ANXX/calls/Sample_NCI0002tumor_E_C4PJ0ANXX.strelka.indels.annotated.txt NCI0002/NCI0002/calls/NCI0002.hapcaller.annotated.txt NCI0002/NCI0002/calls/NCI0002.platypus.annotated.txt NCI0002/NCI0002/calls/NCI0002.bam2mpg.annotated.txt NCI0002/Sample_NCI0002tumor_T_D1UTYACXX/calls/Sample_NCI0002tumor_T_D1UTYACXX.hapcaller.snpEff.txt

import sys
import re
import os 

caller_dict = {}
qual = {}
fs = {}
reads = {}

last_anno_col = 'ACMG_LSDB'

for file_path in sys.argv[1:]:
    a = os.path.basename(file_path).split('.')
    caller = a[1]
    samples = {}
    with open(file_path,'r') as f:
        line_num = 1
        for line in f.readlines():
            line = line.rstrip()
            fields = line.split('\t')
            if line_num == 1 and line[:3] =='Chr':
                header = fields
                last_anno_col_num = header.index(last_anno_col)
                qual_index = header.index('QUAL')
                fs_index = header.index('QUAL')
                for i in range(last_anno_col_num + 5, len(fields), 5):
                    tmp = fields[i]
                    tmp = tmp[:tmp.find('.GT')]
                    samples[i] = tmp
                line_num += 1
                continue   
            if '_' in fields[0]:
                continue
            var = '\t'.join(fields[:5])
            for sample_key in samples.keys():
                if not var + '\t' + samples[sample_key] in reads.keys():
                    reads[var + '\t' + samples[sample_key]] = fields[sample_key +1] + '\t' + fields[sample_key + 3]
                    qual[var + '\t' + samples[sample_key]] = fields[qual_index]
                    fs[var + '\t' + samples[sample_key]] = fields[fs_index]
                    caller_dict[var + '\t' + samples[sample_key]] = caller
                else:
                    qual[var + '\t' + samples[sample_key]] = qual[var + '\t' + samples[sample_key]] +';'+ fields[qual_index]
                    fs[var + '\t' + samples[sample_key]] = fs[var + '\t' + samples[sample_key]] +';'+ fields[fs_index]
                    caller_dict[var + '\t' + samples[sample_key]] = caller_dict[var + '\t' + samples[sample_key]] +';'+ caller
            line_num += 1



for key in sorted(reads.keys()):
    total, alt = reads[key].split('\t')
    if re.search('^\d+$',total) and re.search('^\d+$',alt):
        print('\t'.join([key,caller_dict[key],qual[key],fs[key],reads[key]]))
                                



