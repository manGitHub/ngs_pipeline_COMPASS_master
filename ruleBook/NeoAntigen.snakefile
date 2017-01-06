############
#	seq2HLA
############
rule seq2HLA:
	input:	R=lambda wildcards: FQ[wildcards.sample],
		script=NGS_PIPELINE + "/scripts/seq2HLA.py"
	output:
		"{base}/{TIME}/{sample}/HLA/seq2HLA/{sample}-ClassI.HLAgenotype4digits"
	params:
		rulename= "seq2HLA",
		batch	= config[config['host']]["job_annovar"],
		HLA	= config['HLA']
	shell: """
	#######################
	module load bowtie/1.1.1 python/2.7.10
	python {input.script} {params.HLA}/seq2HLA/ -1 {input.R[0]} -2 {input.R[1]}  -p 2 -r {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/HLA/seq2HLA/ 
	#######################
	"""
############
#	HLAminer
############
rule HLAminer:
	input:	R=lambda wildcards: FQ[wildcards.sample],
		script=NGS_PIPELINE+"/HLAminer_v1.3.1/bin/HPTASRwgs_classI.sh"
	output: "{base}/{TIME}/{sample}/HLA/HLAminer/HLAminer_HPTASR.csv"
	params:
		rulename="HLAminer",
		batch	=config[config['host']]["job_annovar"],
		HLA	=config['HLA'],
		location=NGS_PIPELINE
	shell: """
	#######################
	echo {input.R[0]} >{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/HLA/HLAminer/patient.fof
	echo {input.R[1]} >>{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/HLA/HLAminer/patient.fof
	
	sh {input.script} {params.location}/HLAminer_v1.3.1/bin/ {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/HLA/HLAminer/
		
	#######################
	"""
############
#	pVACSeq
############
rule pVACSeq:
	input:
		files=lambda wildcards: UNION_SOM_MUT[wildcards.sample],
		tool =NGS_PIPELINE + "/scripts/consensusSomaticVCF.pl",
		hla_1="{base}/{TIME}/{sample}/HLA/seq2HLA/{sample}-ClassI.HLAgenotype4digits",
		hla_2="{base}/{TIME}/{sample}/HLA/HLAminer/HLAminer_HPTASR.csv",
		merge=NGS_PIPELINE + "/scripts/consensusHLA.pl"
	output: 
		vcf="{base}/{TIME}/{sample}/NeoAntigen/{sample}.somatic.vep.vcf"
	version: config["R"]
	params:
		rulename = "pVACSeq",
		normal	 = lambda wildcards: config['sample_references'][wildcards.sample][0],
		batch    = config[config['host']]["job_VEP"]
	shell: """
	#######################
	module load vcftools VEP pvacseq python/2.7.10
	perl {input.tool} -vcf {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.strelka.indels.raw.vcf,{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.strelka.snvs.raw.vcf,{wildcards.base}/{wildcards.TIME}/{wildcards.sample}/calls/{wildcards.sample}.MuTect.raw.vcf -order {params.normal},{wildcards.sample} -filter REJECT |vcf-subset -u -c {wildcards.sample} >{output.vcf}.tmp
	variant_effect_predictor.pl -i {output.vcf}.tmp --plugin Downstream --plugin Wildtype --terms SO --offline --cache --dir_cache $VEPCACHEDIR --assembly GRCh37 --output_file {output.vcf} --vcf --force_overwrite
	rm -rf {output.vcf}.tmp
	
	perl {input.merge} {input.hla_2} {input.hla_1}	>{wildcards.base}/{wildcards.TIME}/{params.normal}/HLA/{params.normal}.Calls.txt

	allele=`cut -f1 {wildcards.base}/{wildcards.TIME}/{params.normal}/HLA/{params.normal}.Calls.txt|grep -v Allele|tr '\\n' ',' |sed -e 's/,$//g'`
	pvacseq run -e 8,9,10,11 {output.vcf} {wildcards.sample} ${{allele}} {{NNalign,NetMHC,NetMHCIIpan,NetMHCcons,NetMHCpan,PickPocket,SMM,SMMPMBEC,SMMalign}} {wildcards.base}/{wildcards.TIME}/{wildcards.sample}/NeoAntigen/
	#######################
	"""