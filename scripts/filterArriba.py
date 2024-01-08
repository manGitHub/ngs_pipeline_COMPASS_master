import sys
import pandas as pd

inFile = sys.argv[1]
outFile = sys.argv[2]
fuseCatFile = sys.argv[3]
starFusFile = sys.argv[4]

#breakpoints in fusion catcher and star fusion
fuseCat = pd.read_csv(fuseCatFile,delimiter='\t')
fuseCat = fuseCat["Fusion_point_for_gene_1(5end_fusion_partner)"].values.tolist() + fuseCat["Fusion_point_for_gene_2(3end_fusion_partner)"].values.tolist()
fuseCat = ["chr" + ":".join(x.split(':')[:-1]) for x in fuseCat]
starFus = pd.read_csv(starFusFile,delimiter='\t')
starFus = starFus["RightBreakpoint"].values.tolist() + starFus["LeftBreakpoint"].values.tolist()
starFus = [":".join(x.split(':')[:-1]) for x in starFus]
genes = set(fuseCat + starFus)

# create a column that appends 3 columns of read counts
# filter out fusions with more than 1 zero 
# keep if contains breakpoint in either STARfusion or fusionCatcher
tbl = pd.read_csv(inFile,delimiter='\t')
tbl['allThree'] = tbl[['split_reads1','split_reads2','discordant_mates']].apply(
    lambda x: x.dropna().tolist(),
    axis=1
)
tbl['countZeros'] = [x.count(0) for x in tbl['allThree']]
tbl = tbl.loc[(tbl.countZeros<2) | (tbl['breakpoint1'].isin(genes)) | (tbl['breakpoint2'].isin(genes))]
tbl.to_csv(outFile,index=False,sep='\t')
