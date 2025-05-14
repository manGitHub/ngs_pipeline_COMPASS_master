# usage:
# sigProfiler Assignment script for use in ngs pipeline

import sys
import pandas as pd
import os
from SigProfilerAssignment import Analyzer as Analyze
import matplotlib.pyplot as plt

vcf = sys.argv[1]
libID = sys.argv[2]
pdf_out = sys.argv[3]

# output for sigPro
samples = libID + '_input'
if not os.path.exists(samples):
    os.mkdir(samples)
os.system('cp ' + vcf + ' ' + samples)
# run sigPro
output = libID + '_output/'
Analyze.cosmic_fit(samples, output, input_type="vcf", context_type="96", cosmic_version=3.4, exome=True,genome_build="GRCh37", signature_database=None,exclude_signature_subgroups=None, export_probabilities=True,export_probabilities_per_mutation=True, make_plots=True,sample_reconstruction_plots='both', verbose=True)
os.system('cp ' + output + '/Assignment_Solution/Activities/Assignment_Solution_Activity_Plots.pdf ' + libID + '_Assignment_Solution_Activity_Plots.pdf')
# recreate plot
solut = output + '/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
solut = pd.read_table(solut, sep='\t')
solut = solut.loc[:, (solut != 0).all()]
x = solut.iloc[0][1:].tolist()
y = solut.columns[1:].tolist()
print(x)
print(y)
plt.pie(x, labels = y)
plt.title(libID + ' − Mutational Signatures − v3.4')
plt.show()
plt.savefig(libID + '.sigProfiler.pdf')
# copy pie chart to Actionable dir, remove intermediate dirs
print('Copying pie chart to actionable dir')
os.system('cp ' + libID + '.sigProfiler.pdf ' + pdf_out)
print('Removing ' + output)
print('Removing ' + samples)
os.system('rm -r ' + output)
os.system('rm -r ' + samples)
