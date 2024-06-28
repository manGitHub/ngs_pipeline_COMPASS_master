# To replace GeneAnnotation.pl
# Written by Rob Schultz 9/23
#  Takes in the temporary output from CombineAnnotations.py, Adds ACMG annotations

import sys

#ANNOTATION = hg19_ACMG_v2.txt 
#acmg_path = sys.argv[1]
#acmg_df = pd.read_csv(acmg_path, encoding = 'windows-1252',sep = '\t')
#if 'Gene.refGene' in acmg_df.columns:
#    acmg_df.rename({'Gene.refGene':'ACMG_Gene','Gene_refGene':'ACMG_Gene'}, axis = 1)
#acmg_df['Gene_refGene'].apply(lambda s: s[1:] if s[0] == '#' else s)

#FILE = AnnotationInput.annotations.final.txt.tmp
#tmp_anno_path = sys.argv[2]
#tmp_anno_path = "AnnotationInput.annotations.final.txt.tmp"

#tmp_anno_df = pd.read_csv(tmp_anno_path, sep = '\t')
#out_path = sys.argv[3]

#ANNOTATION = hg19_ACMG_v2.txt 
acmg_path = sys.argv[1]
dict1 = {}
with open(acmg_path,'r',encoding='windows-1252') as a:
    for line in a.readlines():
        if line[-1] == '\n':
            line = line[:-1]
        if line[0] == '#':
            line = line[1:]
        fields = line.split('\t')
        cols = len(fields)
        if 'refGene' in line:
            if 'Gene_refGene' in line:
                line = line.replace('Gene_refGene','ACMG_Gene')
            if 'Gene.refGene' in line:
                line = line.replace('Gene.refGene','ACMG_Gene')
            dict1['Gene_refGene'] = line
            dict1['Gene.refGene'] = line
        else:
            dict1[fields[0]] = line

#FILE = AnnotationInput.annotations.final.txt.tmp
tmp_anno_path = sys.argv[2]
#tmp_anno_path = "AnnotationInput.annotations.final.txt.tmp"
'''
out_path = sys.argv[3]
with open(tmp_anno_path, 'r') as t:
    with open(out_path, 'w') as outfile:
        str1 = '\t-' * cols
        for line in t.readlines():
            if line[-1] == '\n':
                line = line[:-1]
            fields = line.split('\t')
            if len(fields) > 6:
                val1 = fields[6]
                if val1 in dict1.keys():
                    outfile.write(line + '\t' + dict1[val1] + '\n')
                else:
                    outfile.write(line + str1 + '\n')
            else:
                outfile.write(line + str1 + '\n')

'''

with open(tmp_anno_path, 'r') as t:
    #with open(out_path, 'w') as outfile:
        str1 = '\t-' * cols
        for line in t.readlines():
            if line[-1] == '\n':
                line = line[:-1]
            fields = line.split('\t')
            if len(fields) > 6:
                val1 = fields[6]
                if val1 in dict1.keys():
                    print(line + '\t' + dict1[val1])
                else:
                    print(line + str1)
            else:
                print(line + str1)

        
































