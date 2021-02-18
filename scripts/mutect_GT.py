import sys

# removes GT from mutect vcf
# usage:
# python3 mutect_GT.py <sample>.normal.vcf <sample>.normal_noGT.vcf

inFile = sys.argv[1]
outFile = open(sys.argv[2],'w')

header = []
dat = []
with open (inFile) as f:
    for each in f:
        if each.startswith('#'):
            header.append(each)
        else:
            sepEach = each.split('\t')
            sepEach[8] = ':'.join(sepEach[8].split(':')[1:])
            sepEach[9] = ':'.join(sepEach[9].split(':')[1:])
            each = '\t'.join(sepEach)
            dat.append(each)

together = header + dat
for x in together:
    outFile.write(x.strip() + '\n')
