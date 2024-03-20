import sys
import pandas as pd

inFile = sys.argv[1]
outFile = sys.argv[2]
fuseCatFile = sys.argv[3]
starFusFile = sys.argv[4]

#breakpoints in fusion catcher and star fusion
fuseCat = pd.read_csv(fuseCatFile,delimiter='\t')
fuseCat['Fusion_point_for_gene_1(5end_fusion_partner)'] = ['chr' + ':'.join(x) for x in fuseCat["Fusion_point_for_gene_1(5end_fusion_partner)"].str.split(':').str[:-1]]
fuseCat['Fusion_point_for_gene_2(3end_fusion_partner)'] = ['chr' + ':'.join(x) for x in fuseCat["Fusion_point_for_gene_2(3end_fusion_partner)"].str.split(':').str[:-1]]
fuseCat['fusions'] = fuseCat['Fusion_point_for_gene_1(5end_fusion_partner)'] + '-' + fuseCat["Fusion_point_for_gene_2(3end_fusion_partner)"]
fuseCat['reverseFusions'] = fuseCat["Fusion_point_for_gene_2(3end_fusion_partner)"] + '-' + fuseCat['Fusion_point_for_gene_1(5end_fusion_partner)']
fuseCat = fuseCat['fusions'].values.tolist() + fuseCat['reverseFusions'].values.tolist()

starFus = pd.read_csv(starFusFile,delimiter='\t')
starFus['RightBreakpoint'] = [':'.join(x) for x in starFus["RightBreakpoint"].str.split(':').str[:-1]]
starFus['LeftBreakpoint'] = [':'.join(x) for x in starFus["LeftBreakpoint"].str.split(':').str[:-1]]
starFus['fusions'] = starFus['RightBreakpoint'] + '-' + starFus['LeftBreakpoint']
starFus['reverseFusions'] = starFus['LeftBreakpoint'] + '-' + starFus['RightBreakpoint']
starFus = starFus['fusions'].values.tolist() + starFus['reverseFusions'].values.tolist()

fuss = set(fuseCat + starFus)

# create a column that appends 3 columns of read counts
# filter out fusions with more than 1 zero 
# keep if contains fusion in either STARfusion or fusionCatcher
tbl = pd.read_csv(inFile,delimiter='\t')
tbl['allThree'] = tbl[['split_reads1','split_reads2','discordant_mates']].apply(
    lambda x: x.dropna().tolist(),
    axis=1
)
tbl['countZeros'] = [x.count(0) for x in tbl['allThree']]
tbl['fusions'] = tbl['breakpoint1'] + '-' + tbl['breakpoint2']
tbl = tbl.loc[(tbl.countZeros<2) | (tbl['fusions'].isin(fuss))]
tbl.to_csv(outFile,index=False,sep='\t')
