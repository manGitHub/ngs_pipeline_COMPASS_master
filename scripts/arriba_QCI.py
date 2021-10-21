import sys
import pandas as pd

# converts arriba output to QCI format

headerFile = sys.argv[1]        # TSO AllFusions.csv header
inFile = sys.argv[2]            # arriba-fusions.txt
genesFile = sys.argv[3]
outFile = sys.argv[4]
countFile = sys.argv[5]

# TSO AllFusions.csv header and columns
header = []
cols = []
with open(headerFile) as f:
    for each in f:
        if each.startswith('#'):
            header.append(each.strip())
        else:
            cols.append(each.strip())
cols = ''.join(cols).split(',')
tso = pd.DataFrame(columns=cols)

# fill df with arriba-fusion.txt input
arriba = pd.read_csv(inFile,delimiter='\t')
tso['Gene A'] = arriba['#gene1']
tso['Gene B'] = arriba['gene2']
tso['Gene A Breakpoint'] = arriba['breakpoint1']
tso['Gene B Breakpoint'] = arriba['breakpoint2']
tso['Filter'] = 'PASS'
tso['KeepFusion'] = 'TRUE'
tso['Alt Split Dedup'] = arriba['split_reads1'] + arriba['split_reads2']
tso['Alt Pair Dedup'] = arriba['discordant_mates']
tso['Fusion Directionality Known'] = 'TRUE'
tso['Caller'] = 'Arriba'

# filter genes
with open(genesFile) as f:
    genes = [each.strip() for each in f]
tso = tso.loc[tso['Gene A'].isin(genes) | tso['Gene B'].isin(genes)]

# write out count
count = open(countFile,'w')
count.write(str(len(tso)))
count.close()

# write out header and df
qci = open(outFile,'w')
for heads in header:
    qci.write(heads + '\n')
tso.to_csv(qci,mode='a',index=False)
qci.close()
