import sys
import pandas as pd
import os

lohhla_output = sys.argv[1]
sequenza_geneLevel = sys.argv[2]
seq_alternateSol = sys.argv[3]

lohhla = pd.read_csv(lohhla_output,sep='\t')
purFlag = 0
if os.stat(sequenza_geneLevel).st_size == 0:
    purFlag = 1
else:
    sequenza = pd.read_csv(sequenza_geneLevel,sep='\t')
purity = pd.read_csv(seq_alternateSol,sep='\t')

def get_lohhla(gene):
    dat = lohhla[lohhla['class'] == gene]
    maxx = float(dat['HLA_type1copyNum_withBAF'][dat['HLA_type1copyNum_withBAF']==dat['HLA_type1copyNum_withBAF'].max()].to_string(index=False))
    minn = float(dat['HLA_type1copyNum_withBAF'][dat['HLA_type1copyNum_withBAF']==dat['HLA_type1copyNum_withBAF'].min()].to_string(index=False))
    highh = dat['HLA_A_type1'][dat['HLA_type1copyNum_withBAF']==dat['HLA_type1copyNum_withBAF'].max()].to_string(index=False)
    loww = dat['HLA_A_type1'][dat['HLA_type1copyNum_withBAF']==dat['HLA_type1copyNum_withBAF'].min()].to_string(index=False)
    pvall = dat['PVal'].iloc[0]
    cnRatio = maxx/minn if minn > 0 else float("inf")
    return(maxx,minn,highh,loww,pvall,cnRatio)

# if sequenza gene level file is empty, this indicates purity < 0.3 and purFlag == 1
def get_seq(gene):
    if purFlag == 1:
        return('lowPurity')
    else:
        dat = sequenza[sequenza['Gene'] == gene]
        bb = dat['B'].to_string(index=False)
    return(bb)

def getLOHStatus(purity, maxAllele, minAllele, HighAllele, LowAllele, cnB, pValue, CNRatio):
    # I want to exclude samples from the get go if
    # we don't have sequenza data or the purity is below .3
    if cnB is None:
        return ("Inconclusive; No Sequenza Data", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    elif cnB == 'lowPurity':
        return ("Inconclusive; purity Level less than .3", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    elif purity < 0.30:
        return ("Inconclusive; purity Level less than .3", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    elif pValue is None:
        return ("Inconclusive; LOHHLA did not produce Pvalue for allelic difference", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    # since we have seqeunza data and proper purity level
    # we will try to make some calls 1st LOH
    cnB = [int(float(x.strip())) for x in cnB.split(',')]
    if CNRatio >= 2 and maxAllele-minAllele >= .6 and pValue <= 0.05 and min(cnB) == 0:
        return ("LOH", "{} Lost".format(LowAllele), "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    # Now we will try allelic imbalance
    elif CNRatio >= 1.5 and maxAllele-minAllele >= .6 and pValue <= 0.05 and min(cnB) != 0 and maxAllele >=.7 and minAllele >=.7 :
        return ("Allelic Imbalance", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    # If everything hits for LOH but not the Sequenza then report inconclusive results don't agree with sequenza
    elif CNRatio >= 2 and maxAllele-minAllele >= .6 and pValue <= 0.05 and min(cnB) != 0:
        return ("Inconclusive; Results Don't agree with Sequenza", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    # No LOH there seem to be a bunch of
    # situations we'd call no LOH most obvious
    elif CNRatio <= 1.2 and min(cnB) !=0 and pValue > 0.05:
        return ('No LOH', '-', "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    # if its greater than 1.2 but still not rising to allelic imbalance
    elif  1.2 <= CNRatio < 1.5 and min(cnB)!=0 and maxAllele >=.7 and minAllele >=.7  and maxAllele-minAllele <= 1 and pValue >0.05:
        return ('No LOH', '-', "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    # if after all of that Sequenza is reporting a loss but LOHHLA isn't
    # include that in the note
    elif min(cnB)==0:
        return ("Inconclusive; Sequenza shows loss LOHHLA doesn't predict one", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    elif CNRatio >= 1.5 and maxAllele-minAllele >= .6  and min(cnB) != 0 and maxAllele >=.7 and minAllele >=.7 :
        return ("Inconclusive; LOHHLA doesn't show significant Pvalue for copy difference between alleles", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))
    else:
        return ("Inconclusive; did not meet criteria for calling LOH/No LOH or Allelic Imbalance", "-", "{}({}),{}({})".format(LowAllele, minAllele, HighAllele, maxAllele,))

pur = float(purity['cellularity'][purity['SLPP']==purity['SLPP'].max()].to_string(index=False))

# HLA-A
try:
    Amax,Amin,Ahigh,Alow,ApVal,AcnR = get_lohhla('hla_a')
except:
    print('\nLOH status for HLA-A:\nNo LOHHLA results')
aB = get_seq('HLA-A')
try:
    hla_a = getLOHStatus(pur,Amax,Amin,Ahigh,Alow,aB,ApVal,AcnR)
    print('\nLOH status for HLA-A:')
    print(hla_a)
except Exception:
    pass

# HLA-B
try:
    Bmax,Bmin,Bhigh,Blow,BpVal,BcnR = get_lohhla('hla_b')
except:
    print('\nLOH status for HLA-B:\nNo LOHHLA results')
bB = get_seq('HLA-B')
try:
    hla_b = getLOHStatus(pur,Bmax,Bmin,Bhigh,Blow,bB,BpVal,BcnR)
    print('\nLOH status for HLA-B:')
    print(hla_b)
except Exception:
    pass

# HLA-C
try:
    Cmax,Cmin,Chigh,Clow,CpVal,CcnR = get_lohhla('hla_c')
except:
    print('\nLOH status for HLA-C:\nNo LOHHLA results')
cB = get_seq('HLA-C')
try:
    hla_c = getLOHStatus(pur,Cmax,Cmin,Chigh,Clow,cB,CpVal,CcnR)
    print('\nLOH status for HLA-C:')
    print(hla_c)
except Exception:
    pass
