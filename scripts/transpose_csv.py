import pandas as pd
import sys
import os

location = sys.argv[1]
host = sys.argv[2]
version = sys.argv[3]
csv = sys.argv[4]

# if RNAseq, script will take an extra argument with arriba count
# if count == 0, print FUSION NEGATIVE below
count = False
if len(sys.argv) == 6:
	arriba = sys.argv[5]
	with open(arriba) as f:
		count = ''.join([each.strip() for each in f])

# read in csv and reformat 
f = pd.read_csv(csv,header=None,sep='\t')
f = f.T.to_csv(sep=':',index=False,header=False)
f = f.replace('#Patient','Patient').replace(':',': ')

print('Hello,\n')
print('NGS-pipeline version ' + version + ' finished successfully on '+ host + '.\n')
print(f)
if count:
	if count == '0':
		print('FUSION NEGATIVE\n')
print('Results available in ' + location + '.\n')
print('Regards,\nCOMPASS Program\nLaboratory of Pathology\nCCR NCI NIH')
