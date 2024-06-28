# To replace AddSampleType.pl
# Written by Rob Schultz 10/23
# This script creates a temporary output that adds the sample type and capture # method to a vcf-like annotation

import sys

vcf_in = sys.argv[1]
samples = sys.argv[2].split(' ')
sample_dict = {}

for i in range(1, len(samples), 2):
    sample_dict[samples[i-1]] = samples[i]

capture_methods = sys.argv[3].split(' ')
capture_dict = {}
for i in range(1,len(capture_methods),2):
    capture_dict[capture_methods[i-1]] = capture_methods[i]


with open(vcf_in, 'r') as v:
    for line in v.readlines():
        line = line.rstrip()
        fields = line.split('\t')
        print('\t'.join(fields[0:6] + [sample_dict[fields[5]], capture_dict[fields[5]]] + fields[6:11]))
    
