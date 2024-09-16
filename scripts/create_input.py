import pandas as pd
import sys
import os

# usage:
# python3 create_input.py MD# tumorBam normalBam HLA-HD sequenza_alternative_solutions.txt

mrn = sys.argv[1]
tumor_bam_location_pipeline = sys.argv[2]
normal_bam_location_pipeline = sys.argv[3]
HD = sys.argv[4]
seq = sys.argv[5]
LOHHLA_dir = sys.argv[6]

out_file = LOHHLA_dir + '/SnakeMake_input.txt'
tumor_sample_name = os.path.basename(tumor_bam_location_pipeline).split('.bwa.final.bam')[0]
normal_sample_name = os.path.basename(normal_bam_location_pipeline).split('.bwa.final.bam')[0]

CP = tumor_sample_name.split('_')[0]
tumor_bam_location_LOHHLA = LOHHLA_dir+'/'+tumor_sample_name+'_recal.bam'
normal_bam_location_LOHHLA = LOHHLA_dir+'/'+normal_sample_name+'_recal.bam'

# transfer and rename bams
if os.path.isfile(LOHHLA_dir+'/'+tumor_sample_name+'_recal.bam'):
    print('\nTumor BAM was previously transferred.\n')
else:
    print('\nTransferring tumor BAM: '+tumor_sample_name)
    print('cp ' + tumor_bam_location_pipeline + ' ' + tumor_bam_location_LOHHLA)
    print('cp ' + tumor_bam_location_pipeline + '.bai ' + tumor_bam_location_LOHHLA+'.bai')
    print('\n')
    os.system('cp ' + tumor_bam_location_pipeline + ' ' + tumor_bam_location_LOHHLA)
    os.system('cp ' + tumor_bam_location_pipeline + '.bai ' + tumor_bam_location_LOHHLA+'.bai')
if os.path.isfile(LOHHLA_dir+'/'+normal_sample_name+'_recal.bam'):
    print('\nNormal BAM was previously transferred.\n')
else:
    print('\nTransferring normal BAM: '+normal_sample_name)
    print('cp ' + normal_bam_location_pipeline + ' ' + normal_bam_location_LOHHLA)
    print('cp ' + normal_bam_location_pipeline + '.bai ' + normal_bam_location_LOHHLA+'.bai')
    print('\n')
    os.system('cp ' + normal_bam_location_pipeline + ' ' + normal_bam_location_LOHHLA)
    os.system('cp ' + normal_bam_location_pipeline + '.bai ' + normal_bam_location_LOHHLA+'.bai')

# parse HLAs
HD = open(HD,'r')
HD = [line.split('\t') for line in HD]
HD = {line[0] : [item.strip() for item in line[1:]] for line in HD}
HD = pd.DataFrame.from_dict(HD, orient="index")
HD.loc[HD.iloc[:,1]=='-',1] = HD.iloc[:,0]
HD = HD.unstack().to_frame()
HD.columns = ["#Allele"]
HD = HD.dropna()
classI = HD[HD["#Allele"].str.contains("HLA-A|HLA-B|HLA-C",regex=True)]
classI = '|'.join(classI['#Allele'].str.replace('-','_',regex=False).str.replace(':','_',regex=False).str.replace('*','_',regex=False).str.lower().to_list())
classII = HD[HD["#Allele"].str.contains("HLA-DRB1|HLA-DQA1|HLA-DQB1",regex=True)]
classII = '|'.join(classII['#Allele'].str.replace('-','_',regex=False).str.replace(':','_',regex=False).str.replace('*','_',regex=False).str.lower().to_list())

# tumor purity
seq = pd.read_csv(seq,sep='\t')
pur = seq['cellularity'][seq['SLPP']==seq['SLPP'].max()].to_string(index=False)
plo = seq['ploidy'][seq['SLPP']==seq['SLPP'].max()].astype(int)

dat = pd.DataFrame(columns=['mrn','sample_center','Matched_Normal_used','tumor_type','seq center','tumor_sample_name','normal_sample_name','tumor_bam_name','normal_bam_name','tumor_bam_location','normal_nam_location','HLALOH_tumor_name','HLALOH_normal_name','ClassI','ClassII','tumor purity','tumor ploidy','allele check I','allele check II'])
dat['mrn'] = [mrn]
dat['sample_center'] = 'Compass'
dat['Matched_Normal_used'] = normal_sample_name
dat['tumor_type'] = 'tumor_tissue'
dat['seq center'] = 'LP'
dat['tumor_sample_name'] = tumor_sample_name
dat['normal_sample_name'] = normal_sample_name
dat['tumor_bam_name'] = tumor_sample_name+'_recal.bam'
dat['normal_bam_name'] = normal_sample_name+'_recal.bam'
dat['tumor_bam_location'] = tumor_bam_location_LOHHLA
dat['normal_nam_location'] = normal_bam_location_LOHHLA
dat['HLALOH_tumor_name'] = tumor_sample_name+'_recal'
dat['HLALOH_normal_name'] = normal_sample_name+'_recal'
dat['ClassI'] = classI
dat['ClassII'] = classII
dat['tumor purity'] = pur
dat['tumor ploidy'] = plo
dat['allele check I'] = 'Class I input passes'
dat['allele check II'] = 'Class II input passes'
dat.to_csv(out_file,index=False,sep='\t')
