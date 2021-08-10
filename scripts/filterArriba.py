import sys
import pandas as pd

inFile = sys.argv[1]
outFile = sys.argv[2]

# create a column that appends 3 columns of read counts
# filter out fusions with more than 1 zero 
tbl = pd.read_csv(inFile,delimiter='\t')
tbl['allThree'] = tbl[['split_reads1','split_reads2','discordant_mates']].apply(
    lambda x: x.dropna().tolist(),
    axis=1
)
tbl['countZeros'] = [x.count(0) for x in tbl['allThree']]
tbl = tbl.loc[tbl.countZeros<2]
tbl.to_csv(outFile,index=False,sep='\t')
