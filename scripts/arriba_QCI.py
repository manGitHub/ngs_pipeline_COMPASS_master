import sys
import pandas as pd

# converts arriba output to QCI vcf format

inFile = sys.argv[1]            # arriba-fusions.txt
genesFile = sys.argv[2]
outFile = sys.argv[3]
countFile = sys.argv[4]
sampleName = sys.argv[5]

arriba = pd.read_csv(inFile,delimiter='\t')

# count all supporting reads
arriba['reads'] = arriba['split_reads1'].astype(int) + arriba['split_reads2'].astype(int) + arriba['discordant_mates'].astype(int)

# filter genes
with open(genesFile) as f:
    genes = [each.strip() for each in f]
arriba = arriba.loc[arriba['#gene1'].isin(genes) | arriba['gene2'].isin(genes)]
arriba = arriba.loc[arriba['reads'].astype(int) >= 10]

# determine ALT column based on strands
arriba['strand1'] = arriba['strand1(gene/fusion)'].str.split('/').str[0]
arriba['strand2'] = arriba['strand2(gene/fusion)'].str.split('/').str[0]
arriba.loc[(arriba.strand1 == '+') & (arriba.strand2 == '+'), 'ALT'] = 'G[' + arriba['breakpoint2'] + '['
arriba.loc[(arriba.strand1 == '-') & (arriba.strand2 == '+'), 'ALT'] = '[' + arriba['breakpoint2'] + '[G'
arriba.loc[(arriba.strand1 == '+') & (arriba.strand2 == '-'), 'ALT'] = 'G]' + arriba['breakpoint2'] + ']'
arriba.loc[(arriba.strand1 == '-') & (arriba.strand2 == '-'), 'ALT'] = ']' + arriba['breakpoint2'] + ']G'
arriba.loc[(arriba.strand1 == '.') & (arriba.strand2 == '-'), 'ALT'] = 'G]' + arriba['breakpoint2'] + ']'
arriba.loc[(arriba.strand1 == '.') & (arriba.strand2 == '+'), 'ALT'] = 'G[' + arriba['breakpoint2'] + '['
arriba.loc[(arriba.strand1 == '-') & (arriba.strand2 == '.'), 'ALT'] = '[' + arriba['breakpoint2'] + '[G'
arriba.loc[(arriba.strand1 == '+') & (arriba.strand2 == '.'), 'ALT'] = 'G[' + arriba['breakpoint2'] + '['

# if intergenic, take the first gene
arriba['#gene1'] = arriba['#gene1'].str.split(',').str[0]
arriba['gene2'] = arriba['gene2'].str.split(',').str[0]

# vcf
vcf = pd.DataFrame(columns=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sampleName])
vcf['#CHROM'] = arriba['breakpoint1'].str.split(':').str[0]
vcf['POS'] = arriba['breakpoint1'].str.split(':').str[1]
vcf['ID'] = '.'
vcf['REF'] = 'G'
vcf['ALT'] = arriba['ALT']
vcf['QUAL'] = '.'
vcf['FILTER'] = 'PASS'
vcf['INFO'] = 'SVTYPE=Fusion;READ_COUNT=' + arriba['reads'].astype(str) + ';GENE_1=' + arriba['#gene1'].astype(str) + ';GENE_2=' + arriba['gene2'].astype(str) + ';'
vcf['FORMAT'] = 'GT:GQ'
vcf[sampleName] = './.:.'

# write out count
count = open(countFile,'w')
count.write(str(len(arriba)))
count.close()

# write out vcf
vcf.to_csv(outFile,index=False)
