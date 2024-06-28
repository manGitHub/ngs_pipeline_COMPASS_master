# To replace mpileup.pl
# Written by Rob Schultz 10/23
# appends RNAseq info to mutect and strelka annotated txt vcf files

import sys
import subprocess
import re
import os

vcf_path = sys.argv[1]
bam_path = sys.argv[2]
rna_path = sys.argv[3]
#output_path = sys.argv[4]

def MySearch(chr1,start,end):
    for indel in expressed.keys():
        ind = indel.split('\t')
        if not ind[0][:3] == 'Chr':
            ind[1], ind[2] = int(ind[1]), int(ind[2])
            if ind[0] == chr1 and (start in range(ind[1]-10,ind[1]+10) or end in range(ind[2]-10,ind[2]+10)):
                return("match", expressed[indel])
    return('not matching', 'NA')


expressed = {}
with open(rna_path,'r') as r:
    for line in r.readlines():
        line = line.rstrip()
        fields = line.split('\t')
        expressed['\t'.join(fields[:5])] = '\t'.join(fields[9:14])

rna_sample_name = os.path.basename(bam_path).split('.')[0]
'''
with open(output_path, 'w') as out:
    with open(vcf_path, 'r') as v:
        for line in v.readlines():
            line = line.rstrip()
            fields = line.split('\t')
            if line[:3] == 'Chr' or line[0] == '#':
                out.write('\t'.join([line, rna_sample_name+'.GT', 'RNASeq.TotCov', 'RNASeq.RefCov', 'RNASeq.VarCov', 'RNASeq.VAF\n']))
            elif fields[1] == fields[2] and len(fields[3]) == len(fields[4]) and len(fields[3]) < 2:
                rnaseq = subprocess.run('samtools mpileup -d1000000000000 -r '+fields[0]+':'+fields[1]+'-'+fields[1]+' ' + bam_path + ' 2>/dev/null |cut -f 3-5',capture_output = True, text = True,shell = True, check = True).stdout
                rnaseq = rnaseq.rstrip()
                if re.search(r'\d', rnaseq): # Have Coverage 
                    sam_fields = rnaseq.split('\t')    
                    if int(sam_fields[1]) < 1:
                        out.write('\t'.join([line,'NA','0','0','0','0\n']))
                    else:
                        sam_fields[2] = sam_fields[2].upper()
                        bases = ''.join(sorted(sam_fields[2]))
                        ref = bases.count(fields[3])
                        alt = bases.count(fields[4])
                        total_ref = ref + alt
                        if total_ref >= 1:
                            if alt >= 1:
                                vaf = alt/total_ref
                                out.write('\t'.join([line,'NA',str(total_ref),str(ref),str(alt),str(vaf) + '\n']))
                            else:
                                out.write('\t'.join([line,'NA',str(total_ref),str(ref),str(alt),'0\n']))
                        else:
                            out.write('\t'.join([line,'NA',str(0),str(0),str(0),'0\n']))
                else: # no output from mpileup
                    out.write('\t'.join([line,'NA',str(0),str(0),str(0),'0\n']))
            else:
                if '\t'.join(fields[:5]) in expressed.keys():
                    out.write('\t'.join([line, expressed['\t'.join(fields[:5])] + '\n']))
                else:
                    status, info = MySearch(fields[0], fields[1],fields[2])  
                    if status == 'match':
                        out.write(line + '\t' + info + '\n')
                    else:
                        out.write('\t'.join([line,'NA',str(0),str(0),str(0),'0\n']))

'''

with open(vcf_path, 'r') as v:
        for line in v.readlines():
            line = line.rstrip()
            fields = line.split('\t')
            if line[:3] == 'Chr' or line[0] == '#':
                print('\t'.join([line, rna_sample_name+'.GT', 'RNASeq.TotCov', 'RNASeq.RefCov', 'RNASeq.VarCov', 'RNASeq.VAF']))
            elif fields[1] == fields[2] and len(fields[3]) == len(fields[4]) and len(fields[3]) < 2:
                rnaseq = subprocess.run('samtools mpileup -d1000000000000 -r '+fields[0]+':'+fields[1]+'-'+fields[1]+' ' + bam_path + ' 2>/dev/null |cut -f 3-5',capture_output = True, text = True,shell = True, check = True).stdout
                rnaseq = rnaseq.rstrip()
                if re.search(r'\d', rnaseq): # Have Coverage 
                    sam_fields = rnaseq.split('\t')    
                    if int(sam_fields[1]) < 1:
                        print('\t'.join([line,'NA','0','0','0','0']))
                    else:
                        sam_fields[2] = sam_fields[2].upper()
                        bases = ''.join(sorted(sam_fields[2]))
                        ref = bases.count(fields[3])
                        alt = bases.count(fields[4])
                        total_ref = ref + alt
                        if total_ref >= 1:
                            if alt >= 1:
                                vaf = alt/total_ref
                                print('\t'.join([line,'NA',str(total_ref),str(ref),str(alt),str(vaf)]))
                            else:
                                print('\t'.join([line,'NA',str(total_ref),str(ref),str(alt),'0']))
                        else:
                            print('\t'.join([line,'NA',str(0),str(0),str(0),'0']))
                else: # no output from mpileup
                    print('\t'.join([line,'NA',str(0),str(0),str(0),'0']))
            else:
                if '\t'.join(fields[:5]) in expressed.keys():
                    print('\t'.join([line, expressed['\t'.join(fields[:5])]]))
                else:
                    status, info = MySearch(fields[0], fields[1],fields[2])  
                    if status == 'match':
                        print(line + '\t' + info)
                    else:
                        print('\t'.join([line,'NA',str(0),str(0),str(0),'0']))

