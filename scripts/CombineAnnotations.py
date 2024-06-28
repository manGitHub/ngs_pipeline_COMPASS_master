# Adapted from CombineAnnotations.pl
# Written by Rob Schultz 9/23
# Combines annotations, specifically AnnotationInput and AnnotationInput.gene.
# Script is generalizable to any set of annotations given to it. 
#
#

import sys

#print(sys.argv)
infile = sys.argv[1]

with open(infile, 'r') as f:
    paths = f.read().strip().split('\n')



all_sites = {}
AnnotationInput_path = paths[0]
with open(AnnotationInput_path,'r') as A:
    for line in A.readlines():
        if line[-1] == '\n':
            line = line[:-1]
        fields = line.split('\t')
        key = '\t'.join(fields[0:5])
        all_sites[key] = ''

for path in paths:
    with open(path, 'r') as anno_file1:
        anno1 = {}
        #cols = 0
        for line in anno_file1.readlines():
            line = line.strip()
            if line[0] == '#':
                line = line[1:]
            fields = line.split('\t')
            key = '\t'.join(fields[0:5])
            end = len(fields)
            cols = end - 5 
            value = '\t'.join(fields[5:end])
            if not key in anno1.keys():
                anno1[key] = value

    # Done with that annotation
    #Add new annotation (anno1) to the full dict (all_sites)
    for key in all_sites.keys():
        if key in anno1.keys():
            all_sites[key] = all_sites[key] + anno1[key] + '\t'
        else:
            all_sites[key] = all_sites[key] + '0\t'*cols


# Final printing and some manipulations
'''with open(sys.argv[2],'w') as o:
    lines = []
    for key in sorted(all_sites.keys()):
        all_sites[key] = all_sites[key].rstrip()
        lines.append(key + all_sites[key] + '\n')
    o.writelines(lines)
'''
for key in sorted(all_sites.keys()):
    all_sites[key] = all_sites[key].rstrip()
    print(key + all_sites[key])



























