import sys
from sys import argv
import os
import numpy as np
import pandas as pd
import pybedtools
from pybedtools import BedTool
import re

#Reading all the input files:
input_file = sys.argv[1]
bed_file_input= sys.argv[2]
gene_list = sys.argv[3]

#Reading all the temp files:
segment_temp1= sys.argv[4]
py_intersection1= sys.argv[5]
result_temp1=sys.argv[6]
gene_seperated_file=sys.argv[7]
py_intersection2 =sys.argv[8]
final_temp= sys.argv[9]

def gene_agg(data):
                 aggs_res = data.sort_values(by=['#chromosome','start.pos','end.pos']).groupby('Gene',sort= False,as_index=False).agg({'#chromosome': 'first', 'start.pos': 'first','end.pos': 'last','CNt': 'median', 'A': 'median','B':'median'})
                 aggs_res['CNt']= aggs_res['CNt'].round(0).astype(int)
                 aggs_res['A'] = aggs_res['A'].astype(int)
                 aggs_res['B'] = aggs_res['B'].astype(int)
                 return aggs_res


#Removing quotes from the input file
with open(input_file, 'r') as f, open(segment_temp1, 'w') as fo:
    headerline=f.readline()
    headerline = headerline.replace('"', '').replace("'", "")
    fo.write("#" + headerline + "\n")
    for line in f:
        fo.write(line.replace('"', '').replace("'", ""))


#bedtools intersection using pybedtools of the bed_interval file and input file without quotations (segment_temp1)
bed_file = pybedtools.BedTool(bed_file_input)
segment_file = pybedtools.BedTool(segment_temp1)

result_temp= bed_file.intersect(segment_file, wa = True)
result = result_temp.saveas(py_intersection1)

#Replacing all the not found genes and extracting the 1st four columns 
with open(py_intersection1, 'r') as f , open (result_temp1, 'w') as final:
         final.write("Chrom"+ '\t' + "Start.pos" + '\t' + "End.pos" + '\t' + "Gene" + '\n')
         list_random = []
         for line in f:
             line_without=line.replace('___', '\t')
             list_random = line_without.split()
             if not "NOTFOUND" in line_without:
                     final.write("\t".join(list_random[0:4]) + '\n')

#Spliting all the multiple genes value to different rows
df = pd.read_csv(result_temp1, delimiter = '\s+', index_col=False)
df= df.drop('Gene', axis=1).join(df['Gene'].str.split(',', expand=True).stack().reset_index(level=1, drop=True).rename('Gene'))
df.to_csv(gene_seperated_file, header=None, index=None, sep='\t')

#Bedtools interesection with the newgene separated file and the segment file
temp1_file = pybedtools.BedTool(gene_seperated_file)

result_intersection = segment_file.intersect(temp1_file,wb = True)
result = result_intersection.saveas(py_intersection2)

#Final temp file with all the bed and list details
line_final= []
with open(py_intersection2, 'r') as f , open (final_temp, 'w') as fo:
         #fo.write('#chromosome\tstart.pos\tend.pos\tGene\tBf\tN.BAF\tsd.BAF\tdepth.ratio\tN.ratio\tsd.ratio\tCNt\tA\tB\tLPP\tGene\tSource;Source\t#Lists\n')
         for line in f:
             list_final = line.split()
             gene_index= list_final[-1:]
             first_index = list_final[0:3]
             middle_index= list_final[3:13]
             fo.write( "\t".join(first_index) + "\t" +"\t".join(gene_index)+"\t"+ "\t".join(middle_index) + '\n')

#Check if the temporary file is empty
if os.stat(final_temp).st_size != 0:

             #Reading the combined gene list text file and converting into a dictionary
             gene_list = sys.argv[3]
             ann_dict = {}
             ann_array =[]
             with open(gene_list) as ann_read:

                 for line in ann_read:
                     ann_array = line.split()
                     key_gene= ann_array[0]
                     ann_dict[key_gene] = ann_array[1], ann_array[2]

             #Annotation of the genes
             temp_array = []
             OutputFile = open(sys.argv[10], 'w')
             OutputFile.write('#chromosome\tstart.pos\tend.pos\tGene\tBf\tN.BAF\tsd.BAF\tdepth.ratio\tN.ratio\tsd.ratio\tCNt\tA\tB\tLPP\tgene_annotation\tSource;Source\t#Lists\n')

             with open(final_temp, 'r') as temp1_read:
             #next(temp1_read)
                 for line in temp1_read:

                      temp_array = line.split()
                      gene_value = temp_array[3]
           
                      if gene_value in ann_dict.keys():
          
                           list_values = list(ann_dict[gene_value])
                           OutputFile.write("\t".join(temp_array) + '\t' + gene_value + '\t' + "\t".join(list_values) + "\n")
                      else:
                           OutputFile.write("\t".join(temp_array) + '\t' + '-' + '\t'+ '-' + '\t'+ '-' + '\n')
             OutputFile.close()
             #aggregating at the gene level and taking the median for CNV
            
             df_read = pd.read_csv(sys.argv[10],delimiter = '\t')
             df_1 = gene_agg(df_read)
             df_final= df_1[['#chromosome','start.pos','end.pos','Gene','CNt', 'A','B']]
             df_final.to_csv(sys.argv[11],sep='\t',header=True, index= False)
else:
    print('The temp file is empty; coming out of the pipeline')


