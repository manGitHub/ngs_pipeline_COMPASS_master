import os
import subprocess
import argparse
import os
import shutil
import re
from pathlib import Path
import subprocess
from subprocess import Popen, PIPE
import pandas as pd
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('tumorbam')
parser.add_argument('normalbam')
parser.add_argument('seqzpurity')
parser.add_argument('hlas')
parser.add_argument('mrn')
parser.add_argument('outdir')
parser.add_argument('patientID')
parser.add_argument('script')
parser.add_argument('hlafasta')
parser.add_argument('hladat')
parser.add_argument('coordinateFile')
parser.add_argument('tableGenScript')
parser.add_argument('lohTable')
parser.add_argument('lohPlot')
parser.add_argument('lohNotProc')
parser.add_argument('rLibs')
args = parser.parse_args()

tumorbam = args.tumorbam
normalbam = args.normalbam
seqzpurity = args.seqzpurity
hlas = args.hlas
mrn = args.mrn
outdir = args.outdir
patientID = args.patientID
script = args.script
hlafasta = args.hlafasta
hladat = args.hladat
coordinateFile = args.coordinateFile
tableGenScript = args.tableGenScript
lohTable = args.lohTable
lohPlot = args.lohPlot
lohNotProc = args.lohNotProc
rLibs = args.rLibs
print(outdir)

jID = os.environ['SLURM_JOBID']

##make directories
tempdir = '/lscratch/{}'.format(jID)
for x in ['output', 'bamfiles', 'data']:
    if not os.path.exists('{}/{}'.format(tempdir,x)):
        Path('{}/{}'.format(tempdir,x)).mkdir(parents=True, exist_ok=True)

##copy purity and hlas
shutil.copy2(hlas, '{}/data/{}'.format(tempdir, os.path.basename(hlas)))
shutil.copy2(seqzpurity, '{}/data/{}'.format(tempdir, os.path.basename(seqzpurity)))

##copy bams and reindex (there was an issue in another pipeline with the indexes being corupt so I just want to avoid that being a problem)
for x in [tumorbam, normalbam]:
    shutil.copy2(x, '{}/bamfiles/{}'.format(tempdir,os.path.basename(x)))
    #cmd = 'module load singularity ; /data/Compass/local/apps/samtools/1.3.1/bin/samtools index {}/bamfiles/{} {}/bamfiles/{}'.format(tempdir, os.path.basename(x), tempdir,os.path.basename(x).replace('.bam','.bai') )
    cmd = 'module load singularity ; samtools index {}/bamfiles/{} {}/bamfiles/{}'.format(tempdir, os.path.basename(x), tempdir,os.path.basename(x).replace('.bam','.bai') )
    p = subprocess.Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    stdout, stderr = p.communicate()
    p.wait()
    print('########INDEX STDERR########')
    print(stderr.decode('utf-8'))
    print('########INDEX STDOUT########')
    print(stdout.decode('utf-8'))

##create Rscript command
#jellyfish = '/data/Compass/local/apps/jellyfish/2.2.10/bin/jellyfish'

## Edited by SS. 02132023. Added version of novocraft to module load novocraft/3.09.02
# cmd0 = 'module load samtools/1.3.1 ; module load novocraft/3.09.02 ; module load jellyfish/2.2.7 ; '#module load bedtools/2.25.0 ;
#cmd0 = 'module load R/4.2.0 ; module load singularity ; '
cmd1 = 'Rscript {} --patientId {} --outputDir {} --normalBAMfile {} --BAMDir {} --hlaPath {} --coordinateFile {} --rLibraries {} '.format(script, os.path.basename(tumorbam).replace('_recal.bam',''),tempdir+'/output', tempdir+'/bamfiles/{}'.format(os.path.basename(normalbam)),tempdir+'/bamfiles', tempdir+'/data/{}'.format(os.path.basename(hlas)), coordinateFile, rLibs)
cmd2 = '--HLAfastaLoc {} --HLAexonLoc {} --CopyNumLoc {} --mappingStep TRUE --coverageStep TRUE '.format(hlafasta, hladat,tempdir+'/data/{}'.format(os.path.basename(seqzpurity)))
#cmd3 = '--minCoverageFilter 30 --fishingStep TRUE --cleanUp FALSE --gatkDir {} --novoDir {} --bedtools {} --jellyfish {} --samtools {}'.format(gatk, novo, bedT, jellyfish, samtools)
cmd3 = '--minCoverageFilter 30 --fishingStep TRUE --cleanUp FALSE --gatkDir {} --bedtools {} --jellyfish {} --samtools {}'.format('gatk', 'bedtools', 'jellyfish', 'samtools')

#Rcmd = cmd0+cmd1+cmd2+cmd3
Rcmd = cmd1+cmd2+cmd3
print(Rcmd)

c = subprocess.Popen(Rcmd, shell=True, stderr=PIPE, stdout=PIPE)
stdout, stderr = c.communicate()
c.wait()
print('########STDERR########')
print(stderr.decode('utf-8'))
print('########STDOUT########')
print(stdout.decode('utf-8'))

##move output plot
cmd = 'cp {}/output/Figures/{}.HLA.pdf {}'.format(tempdir, os.path.basename(tumorbam).replace('.bam',''), lohPlot)
print(cmd)
c = subprocess.Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
stdout, stderr = c.communicate()
c.wait()
print('########MV STDERR########')
print(stderr.decode('utf-8'))
print('########MV STDOUT########')
print(stdout.decode('utf-8'))

cmd = 'module load R/4.2.0 ; Rscript {} --rdatadir {}/output/ --samplename {} --patientID {}'.format(tableGenScript,tempdir,os.path.basename(tumorbam).replace('_recal.bam','') , mrn )
print(cmd)
c = subprocess.Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
stdout, stderr = c.communicate()
c.wait()
print('########Table Create STDERR########')
print(stderr.decode('utf-8'))
print('########Table Create STDOUT########')
print(stdout.decode('utf-8'))

cmd = 'cp {}/output/{}_hla_cn.txt {}'.format(tempdir, os.path.basename(tumorbam).replace('_recal.bam',''), lohTable)
print(cmd)
c = subprocess.Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
stdout, stderr = c.communicate()
c.wait()
print('########MV table STDERR########')
print(stderr.decode('utf-8'))
print('########MV table STDOUT########')
print(stdout.decode('utf-8'))

if os.path.exists('{}/output/Alleles_not_processed.txt'.format(tempdir)):
    print(cmd)
    cmd = 'mv {}/output/Alleles_not_processed.txt {}'.format(tempdir, lohNotProc)
    c = subprocess.Popen(cmd, shell=True, stderr=PIPE, stdout=PIPE)
    stdout, stderr = c.communicate()
    c.wait()
    print('########MV Alleles not processed STDERR########')
    print(stderr.decode('utf-8'))
    print('########MV Alleles not processed STDOUT########')
    print(stdout.decode('utf-8'))
