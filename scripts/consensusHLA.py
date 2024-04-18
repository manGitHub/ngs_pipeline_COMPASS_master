import pandas as pd
import sys
import re

opti = sys.argv[1]
HD = sys.argv[2]
out = sys.argv[3]

# error handling - check that data was created correctly
def checkData(dat):
    dataCheck = 0
    with open(dat) as f:
        if not 'A1' in f.read() and not 'HLA' in f.read():
            dataCheck = 1
    return dataCheck            
        
# optiType 
optiCheck = checkData(opti)
if optiCheck==0:
    opti = pd.read_csv(opti,sep='\t')
    opti = opti.drop(['Unnamed: 0','Reads','Objective'],axis=1)
    opti = opti.transpose()
    opti = opti.dropna()
    opti.columns = ["#Allele"]
    if opti.empty:
        optiCheck=1
    else:
        opti["#Allele"] = 'HLA-' + opti["#Allele"]
        opti["OptiType"] = opti.index + ":Called"
#    opti['OptiType'] = opti.index
    opti = opti.reset_index(drop=True)

# HLA-HD - has variable num of cols
HDcheck = checkData(HD)
if HDcheck==0:
    HD = open(HD,'r')
    HD = [line.split('\t') for line in HD]
    HD = {line[0] : [item.strip() for item in line[1:]] for line in HD}
    HD = pd.DataFrame.from_dict(HD, orient="index")
    HD.loc[HD.iloc[:,1]=='-',1] = HD.iloc[:,0]
    HD = HD.unstack().to_frame()
    HD.columns = ["#Allele"]
    HD = HD.dropna()
    HD = HD[HD["#Allele"].str.contains("HLA-A|HLA-B|HLA-C",regex=True)]
    HD["#Allele"] = HD["#Allele"].str.split(":").str[:2].str.join(":")
    HD["HLA-HD"] = "Called"
    HD = HD.reset_index(drop=True)

# if one fails, write other to output, otherwise merge and output
if optiCheck==1 and HDcheck==1:
    outF = open(out,'w')
    outF.write('#Allele\tOptiType\tHLA-HD')
elif optiCheck==1 and HDcheck==0:
    HD['OptiType'] = 'NA'
    HD = HD[['#Allele','OptiType','HLA-HD']]
    HD.to_csv(out,sep='\t',index=False)
elif optiCheck==0 and HDcheck==1:
    opti['HLA-HD'] = 'NA'
    opti.to_csv(out,sep='\t',index=False)
else:
    # account for homozygous alleles so merge performs as expected
    opti['#Allele'][opti['#Allele'].duplicated()] = opti['#Allele'] + '_dup'
    HD['#Allele'][HD['#Allele'].duplicated()] = HD['#Allele'] + '_dup'
    consen = pd.merge(opti,HD,on="#Allele",how="outer")
    consen['#Allele'] = consen['#Allele'].str.replace('_dup','',regex=False)
    consen = consen.fillna('NA')
    consen.to_csv(out,sep="\t",index=False)
