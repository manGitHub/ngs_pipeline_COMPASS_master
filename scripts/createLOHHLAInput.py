import pandas as pd
import argparse
import itertools

parser = argparse.ArgumentParser()
#parser.add_argument('SnakeInput')
parser.add_argument('TSampName')
parser.add_argument('alleleList')
parser.add_argument('hlaIn')
parser.add_argument('purityIn')
args = parser.parse_args()

LOH_alleles = pd.read_csv(args.alleleList, sep = '\t')
LOH_alleles = LOH_alleles['ID'].tolist()

##read file select entry
#df = pd.read_csv(args.SnakeInput, sep = '\t')
SnakeInput = 'SnakeMake_input.txt'
df = pd.read_csv(SnakeInput,sep='\t')
df = df.loc[df['tumor_sample_name']==args.TSampName].copy()
print(df)

##gather Class I& II, tumor purity & tumor ploidy
cII = df['ClassII'].item()
cI = df['ClassI'].item()
LOHname = df['HLALOH_tumor_name'].item()
pur = df['tumor purity'].item()
cn = df['tumor ploidy'].item()

if len(df.loc[df['allele check I'] == 'Class I input passes']) > 0:
    cI = sorted(list(set(cI.split('|'))))
else:
    cI = []

if len(df.loc[df['allele check II'] == 'Class II input passes']) > 0:
    cII = set(cII.split('|'))
    cII = sorted(['{}_{}_{}'.format(x.split('_')[0], x.split('_')[1].upper(),'_'.join(x.split('_')[2:])) for x in cII])
else:
    cII = []



##write to output
with open(args.hlaIn.format(args.TSampName), 'w') as f:
    ##LOHHLA does only accept some alleles at different resolutions so 
    ##if they have a lower resolution I'm going to check that and use it
    alleles_for_input = []
    for x in cI+cII:
        # print(x)
        allele = '_'.join(x.split('_')[0:4])
        #try 4 digit
        if allele in LOH_alleles:
            alleles_for_input.append(allele)
        elif '{}N'.format(allele) in LOH_alleles:
            alleles_for_input.append('{}N'.format(allele))
        
        ##search for combination starting at 6
        else:
            a = itertools.product(['01','02','03','04'],['01','02','03','04'])
            for prod in a:
                if '{}_{}'.format(allele, prod[0]) in LOH_alleles:
                    alleles_for_input.append('{}_{}'.format(allele, prod[0]))
                    # print('we used this {}_{}'.format(allele, prod[0]))
                    break
                elif '{}_{}N'.format(allele, prod[0]) in LOH_alleles:
                    alleles_for_input.append('{}_{}N'.format(allele, prod[0]))
                    # print('we used this {}_{}'.format(allele, prod[0]))
                    break
                elif '{}_{}_{}'.format(allele, prod[0],prod[1]) in LOH_alleles :
                    alleles_for_input.append('{}_{}_{}'.format(allele, prod[0],prod[1]))
                    # print('we used this {}_{}_{}'.format(allele, prod[0],prod[1]))
                    break
                elif '{}_{}_{}N'.format(allele, prod[0],prod[1]) in LOH_alleles:
                    alleles_for_input.append('{}_{}_{}N'.format(allele, prod[0],prod[1]))
                    # print('we used this {}_{}_{}'.format(allele, prod[0],prod[1]))
                    break
    
    alleles_for_input = '|'.join(alleles_for_input)
    f.write('\n'.join(sorted(alleles_for_input.split('|'))))
    
##Write purity/ploidy to output
pudf = pd.DataFrame(columns=['samp','Ploidy','tumorPurity','tumorPloidy'], data = [[LOHname, 2,pur,cn]])
pudf = pudf.set_index('samp', drop = True)
pudf.to_csv(args.purityIn.format(args.TSampName), sep='\t', index_label=False)
