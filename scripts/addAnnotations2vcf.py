# To replace addAnnotations2vcf.pl
# Written by Rob Schultz 10/23
#  Takes in a completed annotation file i.e. Annotations.final.txt or
#  annotations.coding.rare.txt and attaches the annotations there to the
#  snpEff.txt o

import sys

anno_path = sys.argv[1]
vcf_path = sys.argv[2]
#output_path = sys.argv[3]

anno_dict = {}
# ANN_FH
with open(anno_path, 'r', encoding='windows-1252') as a:
    for line in a.readlines():
        line = line.rstrip()
        fields = line.split('\t')
        key = '\t'.join(fields[:5])
        #end = len(fields) - 1
        value = '\t'.join(fields[5:])
        if not key in anno_dict.keys():
            anno_dict[key] = value

with open(vcf_path, 'r') as v:
    #with open(output_path, 'w') as out:
        for line in v.readlines():
            line = line.rstrip()
            fields = line.split('\t')
            val = '\t'.join(fields[:5])
            vcf = '\t'.join(fields[5:])
            if val in anno_dict.keys():
                print(val + '\t' + anno_dict[val] + '\t' + vcf)
                #out.write(val + '\t' + anno_dict[val] + '\t' + vcf + '\n')





