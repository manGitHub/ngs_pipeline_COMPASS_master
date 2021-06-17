# combine mutect snp and strelka indel tmb

import sys

mutectFile = sys.argv[1]
strelkaFile = sys.argv[2]

with open (mutectFile) as m, open(strelkaFile) as s:
	mut = m.readlines()
	strel = s.readlines()
	
	total_bp = mut[0].split('\t')[1].strip()
	
	non = int(mut[3].split('\t')[1].strip()) + int(strel[3].split('\t')[1].strip())	
	all = int(mut[7].split('\t')[1].strip()) + int(strel[7].split('\t')[1].strip())

	print("Total bases\t" + str(total_bp) + '\n')

	print('Nonsynonymous somatic calls:')
	print('Mutation burden\t' + str(non))
	print('Mutation burden per megabase\t' + '{:.3f}'.format((float(non)/float(total_bp))*1000000)+'\n')

	print('All somatic calls:')
	print('Mutation burden\t' + str(all))
	print('Mutation burden per megabase\t' + '{:.3f}'.format((float(all)/float(total_bp))*1000000)+'\n')
	
