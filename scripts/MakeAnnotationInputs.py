# To replace MakeAnnotationInputs.pl
# Written by Rob Schultz 9/23
# Creates the input files for Annovar gene Annotations from .txt formatted vcf
# Does not create annotation_inputs for pph or sift

import sys

vcf_path = sys.argv[1] 
anno_out_path = vcf_path + '.anno'
#pph_out_path = vcf_path + '.pph'
#sift_out_path = vcf_path + '.sift'

with open(vcf_path,'r') as v:
    with open(anno_out_path,'w') as anno:
#        with open(pph_out_path,'w') as pph:
#            with open(sift_out_path,'w') as sift:
                for line in v.readlines():
                    if line[0] == '#' or line[0:3] == 'Chr' or line[0:3] == 'CHR' or line[0:4] == 'chr\t':
                        anno.write("Chr\tStart\tEnd\tRef\tAlt\n") 
                    else:
                        line = line.strip()
                        fields = line.split('\t')
                        if len(fields) >= 5:
                            anno.write('\t'.join(fields[0:5]) + '\n')
#                        if len(fields[3]) < 2 and len(fields[4]) < 2 and any([i in fields[3] for i in 'ATCG']) and any([i in fields[4] for i in 'ATCG']) and  not '-' in fields[3] and not '-' in fields[4]:
#                            pph.write(fields[0] +':'+fields[1] +'\t'+ fields[3] +'/'+ fields[4] +'\n') 
#                            fields[0] = fields[0].replace('chr','')
#                            sift.write(fields[0] +','+ fields[1] +',1,'+ fields[3] +'/'+ fields[4] +'\n') 















