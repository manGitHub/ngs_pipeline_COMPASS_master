#!/usr/bin/python

import os
import pysam
import re
import string
import sys
import gzip


from numpy import array
from sys import argv
import pandas as pd

from bs4 import BeautifulSoup
import lxml

#Usage:
#python3 combined_python.py mosdepth_thresold_bed qualimal.html mosdepth_region_file hsmetrics_file result_table Patient_id library_id diagnosis
#python3 combined_python.py CP06297_T1D_E_final.thresholds.bed.gz CP06297_T1D_E_qualimap/qualimapReport.html CP06297_T1D_E_final.regions.bed.gz CP06297_N1D_E_format.hsmetrics output_table "CP06297" "CP06297_T1D_E" "N/A"

mosdepth_threshold_file = sys.argv[1]
qualimap_html = sys.argv[2]
mosdepth_region_bed = sys.argv[3]
hsmetrics = sys.argv[4]
final_result = sys.argv[5]




#Working on mosdepth results to get the avearge thereshold values:
stats = {}
def coverage(txt):

    with open (txt, 'rb') as fd:

        gzip_fd = gzip.GzipFile(fileobj=fd)
        df = pd.read_csv(gzip_fd, sep='\t')

        total_positions = 36583011
        stats['total_5X'] = (float(df.iloc[:,4:5].sum()))/ (total_positions) * 100
        stats['total_10X']= (float(df.iloc[:,5:6].sum()))/(total_positions) * 100
        stats['total_15X']= (float(df.iloc[:,6:7].sum()))/(total_positions) * 100
        stats['total_20X']= (float(df.iloc[:,7:8].sum()))/(total_positions) *100
        stats['total_30X']= (float(df.iloc[:,8:9].sum()))/(total_positions) * 100
        stats['total_50X']=(float(df.iloc[:,9:10].sum()))/(total_positions) * 100
        stats['total_100X']=(float(df.iloc[:,10:11].sum()))/(total_positions) * 100
        stats['total_200X']=(float(df.iloc[:,11:12].sum()))/(total_positions) * 100
        stats['total_400X']=(float(df.iloc[:,12:13].sum()))/(total_positions) * 100
        stats['total_1000X']=(float(df.iloc[:,13:14].sum()))/(total_positions) * 100

        return stats
      
def avg_coverage(region_mosdepth):

    with open (region_mosdepth, 'rb') as bed_file:

        gzip_fd2 = gzip.GzipFile(fileobj=bed_file)
        df2 = pd.read_csv(gzip_fd2, delim_whitespace=True)
        stats['average']= float(df2.iloc[:,4:5].mean())         

        return stats['average']



#Working on Qualimap Results:
#html parsing

def parse_html_table(table):
            n_columns = 0
            n_rows=0
            column_names = []
    
            # Find number of rows and columns
            # we also find the column titles if we can
            for row in table.find_all('tr'):
                
                # Determine the number of rows in the table
                td_tags = row.find_all('td')
                if len(td_tags) > 0:
                    n_rows+=1
                    if n_columns == 0:
                        # Set the number of columns for our table
                        n_columns = len(td_tags)
                        
                # Handle column names if we find them
                th_tags = row.find_all('th') 
                if len(th_tags) > 0 and len(column_names) == 0:
                    for th in th_tags:
                        column_names.append(th.get_text())
    
            # Safeguard on Column Titles
            if len(column_names) > 0 and len(column_names) != n_columns:
                raise Exception("Column titles do not match the number of columns")
    
            columns = column_names if len(column_names) > 0 else range(0,n_columns)
            df = pd.DataFrame(columns = columns,
                              index= range(0,n_rows))
            row_marker = 0
            for row in table.find_all('tr'):
                column_marker = 0
                columns = row.find_all('td')
                for column in columns:
                    df.iat[row_marker,column_marker] = column.get_text()
                    column_marker += 1
                if len(columns) > 0:
                    row_marker += 1
            return df 


#Creating the beautiful soup object for parsing:
bs_object = BeautifulSoup(open(qualimap_html, encoding='utf-8'),"lxml")

#Getting tables from the html
table_summary1= bs_object.find_all('table')[3] #for total summary
table_summary2= bs_object.find_all('table')[4] #for target summart
table_mmq = bs_object.find_all('table')[7] #for mean mmq

new_table= parse_html_table(table_summary1)
new_table2= parse_html_table(table_summary2)  
new_table_mmq = parse_html_table(table_mmq)

combined = pd.concat([new_table, new_table_mmq], axis=0, join = "outer", ignore_index = False)
df1 = combined.transpose()
df1.columns= df1.iloc[0]
df1= df1[1:]
df2 = new_table2.transpose()
df2.columns= df2.iloc[0]
df2= df2[1:]


#picard Hsmetrics coverage details
tbl= []
with open (argv[4], 'r') as hs_matrix:
     for each in hs_matrix:
        if not each.startswith('#'):
            tbl.append(each.rstrip())
df_hsmetrics = pd.DataFrame([x.split('\t') for x in tbl])
df_hsmetrics.columns = df_hsmetrics.iloc[0]
df_hsmetrics = df_hsmetrics[1:]
df1_hsmetrics= df_hsmetrics.iloc[:2].reset_index(drop=True)
df1_hsmetrics.columns= df1_hsmetrics.iloc[0]
df1_hsmetrics= df1_hsmetrics[1:]
mean_bait= df1_hsmetrics['MEAN_BAIT_COVERAGE'].item()
mean_target=df1_hsmetrics['MEAN_TARGET_COVERAGE'].item()
print(mean_bait)
print(mean_target)


with open (final_result, 'w') as output:
          
         details= coverage(mosdepth_threshold_file)
         avg_depth=avg_coverage(mosdepth_region_bed)
    
         stats['average']=  "%.2f"% stats['average']       

         Number_of_reads=df1['Number of reads'].item()
         Number_of_reads= Number_of_reads.replace(',' , '')
         Mapped_reads_list= df1['Mapped reads'].str.split(' /')[1]
         Mapped_reads= int(Mapped_reads_list[0].replace(',' , ''))
         Percentage_mapped_reads= Mapped_reads_list[1]

         On_target_reads_list=df2['Mapped reads'].str.split(' /')[1]
         On_target_reads = On_target_reads_list[0].replace(',' , '')
         On_target_percentage = (On_target_reads_list[1])

         Duplicate_read_list = df2['Duplicated reads (flagged)'].str.split(' /')[1] 
         Unique_on_target_reads = int(On_target_reads_list[0].replace(',' , '')) - int(Duplicate_read_list[0].replace(',' , ''))
         Percentage_unique_reads = "%.2f"% ((int(Unique_on_target_reads) / int(On_target_reads)) * 100)         
         mean_mapping_quality = df1['Mean Mapping Quality'].item()

         output.write("Patient\tSample\tDiagnosis\tAverage_coverage_mosdepth\tMEAN_BAIT_COVERAGE\tMEAN_TARGET_COVERGAE\ttotal_reads\tmapped_reads\tpercent_mapped\tontarget_reads\tpercent_ontarget\tunique_on_target_reads\tpercent_unique_on_target_reads\tmean_mapq\tpercent_unique_positions_at_5x\tpercent_unique_positions_at_10x\tpercent_unique_positions_at_15x\tpercent_unique_positions_at_20x\tpercent_unique_positions_at_30x\tpercent_unique_positions_at_50x\tpercent_unique_positions_at_100x\tpercent_unique_positions_at_200x\tpercent_unique_positions_at_400x\tpercent_unique_positions_at_1000x\n")
         output.write((argv[6]) +  '\t' + (argv[7]) + '\t' + (argv[8]) + '\t' + str(stats['average']) + '\t' + str(mean_bait) +'\t' + str(mean_target) + '\t' + str(Number_of_reads) + '\t' + str(Mapped_reads) + '\t' + str(Percentage_mapped_reads) + '\t'  + str(On_target_reads) + '\t' + str(On_target_percentage) + '\t' + str(Unique_on_target_reads) + '\t' + str(Percentage_unique_reads) + '\t' +  str(mean_mapping_quality)+ '\t' + str(stats['total_5X']) + '\t' + str(stats['total_10X']) + '\t' + str(stats['total_15X']) + '\t' + str(stats['total_20X']) + '\t' + str(stats['total_30X']) + '\t' + str(stats['total_50X']) + '\t' + str(stats['total_100X']) + '\t' + str(stats['total_200X']) + '\t' + str(stats['total_400X']) + '\t' + str(stats['total_1000X']))

