import sys

# removes GT from normal mutect vcf after vcfsubset
# usage:
# python3 mutect_GT.py <sample>.normal.vcf <sample>.normal_noGT.vcf

inFile = sys.argv[1]

header = []
dat = []
with open (inFile) as f:
    exitFlag = False
    for each in f:
        if each.startswith('#'):
            header.append(each)
        else:
            sepEach = each.split('\t')
            if sepEach[8].startswith('GT'):
                sepEach[8] = ':'.join(sepEach[8].split(':')[1:])
                sepEach[9] = ':'.join(sepEach[9].split(':')[1:])
                each = '\t'.join(sepEach)
                dat.append(each)
            else:
                exitFlag = True

if not exitFlag:
    outFile = open(sys.argv[2],'w')
    together = header + dat
    for x in together:
        outFile.write(x.strip() + '\n')
else:
    sys.exit('Error: FORMAT field does not begin with GT')
