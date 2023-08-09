import pandas as pd
import sys

opti = sys.argv[1]
HD = sys.argv[2]
out = sys.argv[3]

# read input - HLA-HD has variable num of cols
opti = pd.read_csv(opti,sep='\t')
HD = open(HD,'r')
HD = [line.split('\t') for line in HD]
HD = {line[0] : [item.strip() for item in line[1:]] for line in HD}
HD = pd.DataFrame.from_dict(HD, orient="index")

# convert optitype data
opti = opti.drop(['Unnamed: 0','Reads','Objective'],axis=1)
opti = opti.transpose()
opti.columns = ["#Allele"]
opti["#Allele"] = 'HLA-' + opti["#Allele"]
opti["OptiType"] = "Called"

# convert HLA-HD data
HD = HD.unstack().to_frame()
HD.columns = ["#Allele"]
HD = HD.dropna()
HD = HD[HD["#Allele"].str.contains("HLA-A|HLA-B|HLA-C",regex=True)]
HD["#Allele"] = HD["#Allele"].str.split(":").str[:2].str.join(":")
HD["HLA-HD"] = "Called"

consen = pd.merge(opti,HD,on="#Allele",how="outer")
consen.to_csv(out,sep="\t",index=False)
