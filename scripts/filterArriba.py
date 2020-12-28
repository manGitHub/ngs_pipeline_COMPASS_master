import sys
import pandas as pd

inFile = sys.argv[1]
outFile = sys.argv[2]

# create a column that appends 3 columns of read counts
# filter out fusions with more than 1 zero and low confidence
tbl = pd.read_csv(inFile,delimiter='\t')
tbl['allThree'] = tbl[['split_reads1','split_reads2','discordant_mates']].apply(
    lambda x: ' '.join(x.dropna().astype(str)),
    axis=1
)
tbl['countZeros'] = tbl.allThree.str.count('0')
tbl = tbl.loc[(tbl.countZeros<2) | (tbl.confidence != 'low')]
tbl.to_csv(outFile,index=False,sep='\t')
