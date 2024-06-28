# To replace addFS.pl
# Written by Rob Schultz 10/23
# prints a vcf-like output that contains no header, and each row of the
# input with proper fs information

import sys
import re

vcf_path = sys.argv[1]
FS_dict = {}
for path in sys.argv[2:]:
    if 'HC_DNASeq' in path:
        with open(path,'r') as p:
            for line in p.readlines():
                line = line.rstrip()
                fields = line.split('\t')
                key = '\t'.join(fields[:5] + [fields[8]])
                match = re.search(r'FS=(.*?);',fields[7])
                if match:
                    FS_dict[key] = match.group(1)

            
with open(vcf_path,'r') as v:
    for line in v.readlines():
        line = line.rstrip()
        fields = line.split('\t')
        if 'HC_DNASeq' in fields[8]:
            key = '\t'.join(fields[:6])
            fields[10] = FS_dict[key]
        else:
            fields[10] = 'NA'
        print('\t'.join(fields))

