import sys
import os


ann_var = sys.argv[1]
intervals = sys.argv[2]
tcov = sys.argv[3]
ncov = sys.argv[4]
vaf = sys.argv[5]

# initialize column indices
INDEX_NCOV=200
INDEX_TCOV=205
INDEX_VAF=208
ExonicFunc=8

fname = os.path.basename(ann_var).split('.')[0]

#Calculate total base pairs covered by bed file.
total_bp = 0
with open(intervals) as f:
    iv = f.readlines()
    for each in iv:
        total_bp += int(each.split('\t')[2]) - int(each.split('\t')[1])
    print("Total bases\t" + str(total_bp) + '\n')
f.close()

# filter and count variants
with open(ann_var) as f:
    annvar = f.readlines()
    annvar = annvar[1:] # remove header
    count = 0
    nonsyn = 0
    for each in annvar:
        each = each.split('\t')
        if (float(each[INDEX_NCOV])>=float(ncov)) and (float(each[INDEX_TCOV])>=float(tcov)) and (float(each[INDEX_VAF])>=float(vaf)) and (each[ExonicFunc]not in['-1','unknown']):
            count+=1
            if each[ExonicFunc]not in['synonymous SNV']:
                nonsyn+=1
print('Nonsynonymous somatic calls:')
print('Mutation burden\t' + str(nonsyn))
print('Mutation burden per megabase\t' + '{:.3f}'.format((float(nonsyn)/float(total_bp))*1000000)+'\n')

print('All somatic calls:')
print('Mutation burden\t' + str(count))
print('Mutation burden per megabase\t' + '{:.3f}'.format((float(count)/float(total_bp))*1000000)+'\n')
