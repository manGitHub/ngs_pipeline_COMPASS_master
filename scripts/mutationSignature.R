#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
                make_option("--input", help="input file name"),
                make_option("--sample", help="Header of the sample column"),
		make_option("--output", help="output pdf file name"),
		make_option("--rLib", help="R library location"),
		make_option("--mutSigs", help="mutational signature RData object")
)

opt <- parse_args(OptionParser(option_list=option_list))
input=opt$input
output=opt$output
sample=opt$sample
rLib=opt$rLib
mutSigs=opt$mutSigs

.libPaths(c(rLib, .libPaths()))
suppressPackageStartupMessages(library("deconstructSigs"))

load(mutSigs)

mut_data <- read.table(input,sep="\t",header=T)

sigs.input <- mut.to.sigs.input(mut.ref = mut_data,
		sample.id = "Sample",
		chr = "Chr",
		pos = "Start",
		ref = "Ref",
		alt = "Alt")
# cosmic = whichSignatures(tumor.ref = sigs.input,
#		signatures.ref = signatures.cosmic,
#		sample.id = sample,
#		contexts.needed = TRUE,
#		tri.counts.method = 'default')

# nature = whichSignatures(tumor.ref = sigs.input,
#                signatures.ref = signatures.nature2013,
#               sample.id = sample,
#                contexts.needed = TRUE,
#                tri.counts.method = 'default')

# sig_2019 = whichSignatures(tumor.ref = sigs.input,
#		signatures.ref = signatures.exome.cosmic.v3.may2019,
#		sample.id = sample,
#		contexts.needed = TRUE,
#		tri.counts.method = 'default')

sig_2020 = whichSignatures(tumor.ref = sigs.input,
		signatures.ref = signatures.exome.cosmic.v3.1.june2020_kv,
		sample.id = sample,
		contexts.needed = TRUE,
		tri.counts.method = 'default')

pdf(output,width=10)
#plotSignatures(cosmic, sub='Mutational Signature Based on COSMIC')
#makePie(cosmic, sub='Mutational Signature Based on COSMIC')
#plotSignatures(nature, sub='Mutational Signature based on Nature 2013--23945592')
#makePie(nature, sub='Mutational Signature based on Nature 2013--23945592')
#plotSignatures(sig_2019, sub='Mutational Signatures - 2019 v3')
#makePie(sig_2019, sub='Mutational Signatures - 2019 v3',v3=TRUE)
plotSignatures(sig_2020, sub='Mutational Signatures - 2020 v3.1')
makePie_2020(sig_2020, sub='Mutational Signatures - 2020 v3.1',v3.1=TRUE)

dev.off()
