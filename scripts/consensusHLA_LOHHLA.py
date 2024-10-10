import pandas as pd
import sys
import re

opti = sys.argv[1]
HD = sys.argv[2]
LOHHLA = sys.argv[3]
out = sys.argv[4]

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

# LOHHLA
no_result = ['Inconclusive; No Sequenza Data','Inconclusive; LOHHLA did not produce Pvalue for allelic difference','Allelic Imbalance','Inconclusive; Results Dont agree with Sequenza','Inconclusive; Sequenza shows loss LOHHLA doesnt predict one','Inconclusive; did not meet criteria for calling LOH/No LOH or Allelic Imbalance','No LOHHLA results']
low_purity = ['Inconclusive; purity Level less than .3']
no_LOH = ['No LOH','Inconclusive; LOHHLA doesnt show significant Pvalue for copy difference between alleles']
f = open(LOHHLA,'r')
lohhla = f.readlines()
lohhla = pd.DataFrame([lohhla[2].split(','),lohhla[5].split(','),lohhla[8].split(',')])
if len(lohhla.columns) == 1:
    lohhla['Lost'], lohhla['Allele_1'], lohhla['Allele_2'] = None, None, None
lohhla.columns = ['Status','Lost','Allele_1','Allele_2']
lohhla = lohhla.replace('"','',regex=True).replace("'","",regex=True).replace('\\n','',regex=True)
lohhla['Status'] = lohhla['Status'].str.replace(")","").str.replace("(","")
one = lohhla[['Allele_1','Status','Lost']]
two = lohhla[['Allele_2','Status','Lost']]
one.columns = ['#Allele','Status','Lost']
two.columns = ['#Allele','Status','Lost']
lohhla = pd.concat([one,two],ignore_index=True)
lohhla['LOHHLA'] = ''
lohhla['#Allele'] = lohhla['#Allele'].str.replace("'","").str.split('_').str[0:4].str.join('_').str.replace('_','-',n=1).str.replace('_','*',n=1).str.replace('_',':',n=1).str.split('(').str[0].str.replace(' ','').str.upper()
lohhla['Lost'] = lohhla['Lost'].str.replace("'","").str.split('_').str[0:4].str.join('_').str.replace('_','-',n=1).str.replace('_','*',n=1).str.replace('_',':',n=1).str.split(' Lost').str[0].str.replace(' ','').str.upper()
lohhla['LOHHLA'].loc[lohhla['Lost']!='-'] = lohhla['Lost'] + ' Lost'
lohhla['LOHHLA'].loc[lohhla['Status'].isin(no_result)] = 'No result'
lohhla['LOHHLA'].loc[lohhla['Status'].isin(low_purity)] = 'Low purity'
lohhla['LOHHLA'].loc[lohhla['Status'].isin(no_LOH)] = 'No LOH'
lohhla = lohhla[['#Allele','LOHHLA']]
lohhla = lohhla.dropna()

# if one fails, write other to output, otherwise merge and output
if optiCheck==1 and HDcheck==1:
    outF = open(out,'w')
    outF.write('#Allele\tOptiType\tHLA-HD\tLOHHLA')
elif optiCheck==1 and HDcheck==0:
    HD['OptiType'] = 'NA'
    HD = HD[['#Allele','OptiType','HLA-HD']]
    HD = pd.merge(HD,lohhla,on='#Allele',how='outer')
    HD.to_csv(out,sep='\t',index=False)
elif optiCheck==0 and HDcheck==1:
    opti['HLA-HD'] = 'NA'
    opti['LOHHLA'] = 'NA'
    opti.to_csv(out,sep='\t',index=False)
else:
    # account for homozygous alleles so merge performs as expected
    opti['#Allele'][opti['#Allele'].duplicated()] = opti['#Allele'] + '_dup'
    HD['#Allele'][HD['#Allele'].duplicated()] = HD['#Allele'] + '_dup'
    consen = pd.merge(opti,HD,on="#Allele",how="outer")
    consen['#Allele'] = consen['#Allele'].str.replace('_dup','',regex=False)
    consen = pd.merge(consen,lohhla,on='#Allele',how='outer')
    consen = consen.fillna('NA')
    consen.to_csv(out,sep="\t",index=False)
