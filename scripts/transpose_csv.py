import pandas as pd
import sys

location = sys.argv[1]
host = sys.argv[2]
version = sys.argv[3]
csv = sys.argv[4]
 
f = pd.read_csv(csv,header=None,sep='\t')
f = f.T.to_csv(sep=':',index=False,header=False)
f = f.replace('#Patient','Patient').replace(':',': ')

print('Hello,\n')
print('NGS-pipeline version ' + version + ' finished successfully on '+ host + '.\n')
print(f)
print('Results available in ' + location + '.\n')
print('Regards,\nCOMPASS Program\nLaboratory of Pathology\nCCR NCI NIH')
