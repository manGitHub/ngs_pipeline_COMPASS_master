# To replace UnionSomaticCalls.pl
# Written by Rob Schultz 10/23
# Compiles annotatedFull.txt from Strelka indels and snvs, and MuTect

import sys

#output_path = sys.argv[1]
calls = {}
for path in sys.argv[1:]:
    with open(path,'r') as f: 
        caller = path.split('.')[1]
        for line in f.readlines():
            line = line.rstrip()
            fields = line.split('\t')
            key = '\t'.join(fields[:5])
            if line[:3] == 'Chr':
                calls[key] = line + '\tCaller'
                continue
            calls[key] = line + '\t' + caller


'''
with open(output_path,'w') as out:
    for key in sorted(calls.keys()):
        out.write(calls[key] + '\n')
'''

for key in sorted(calls.keys()):
    print(calls[key])







