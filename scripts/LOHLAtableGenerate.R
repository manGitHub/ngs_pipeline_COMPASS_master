library(optparse)
library(dplyr)

option_list = list(
  make_option(c("-r","--rdatadir"), type="character", default=NULL, 
              help="rdatadir with files to be loaded", metavar="character"),
  make_option(c("-s","--samplename"), type="character", default=NULL, 
              help="name of file header for save", metavar="character"),
  make_option(c("-p","--patientID"), type="character", default=NULL, 
              help="mrn of patient", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


if (is.null(opt$rdatadir)|is.null(opt$samplename)|is.null(opt$patientID)){
  print_help(opt_parser)
  stop("Missing arguments.\n", call.=FALSE)  
}
 
#  print(opt)
rdatadir <- opt$rdatadir
sname <- opt$samplename
patientID <- opt$patientID

df <- data.frame(matrix(ncol=27, nrow=0))
colnames(df) <- c('mrn','sample', 'class','allele', 'Expected_Allele_frequency', 'Observed_allele_frequency', 'LogR_Tumor_vs_Normal' ,
'copy number', 'copyNum_withBAF_lower','copyNum_withBAF_upper',
'copyNum_withBAFBin','copyNum_withBAFBin_lower','copyNum_withBAFBin_upper',
'copyNum_withoutBAF','copyNum_withoutBAF_lower','copyNum_withoutBAF_upper',
'copyNum_withoutBAFBin','copyNum_withoutBAFBin_lower','copyNum_withoutBAFBin_upper', 
'PVal', 'PVal_unique',
'minCoverageFilter','mismatch count','Proportion_of_SupportiveSites',
'total_Mismatch_bases_between_alleles','total_GAPs_between_alleles','Pass_HLA_Allele_Diff_pos_Normal')
print(colnames(df) )

allele_1_mismatch <- NA
allele_1_gaps     <- NA
allele_2_mismatch <- NA 
allele_2_gaps     <- NA


hla_allele_data <- read.table(paste0(rdatadir,"/",sname,".DNA.HLAlossPrediction_CI.txt"), sep="\t",header = TRUE)
if(nrow(hla_allele_data)>1){

all <- do.call(rbind, unlist(apply(hla_allele_data, 1, function(x){
  
  if (!is.na(x["Total_Mismatch_Pos_between_alleles"])) {
    allele_1_mismatch <- unlist(strsplit(x["Total_Mismatch_Pos_between_alleles"],split = " ; "))[1]
    allele_2_mismatch <- unlist(strsplit(x["Total_Mismatch_Pos_between_alleles"],split = " ; "))[2]
  }
  
  if (!is.na(x["Total_GAP_Pos_between_alleles"])) {
    allele_1_gaps <- unlist(strsplit(x["Total_GAP_Pos_between_alleles"],split = " ; "))[1]
    allele_2_gaps <- unlist(strsplit(x["Total_GAP_Pos_between_alleles"],split = " ; "))[2]
  }
  
  hla = paste(unlist(strsplit( x["HLA_A_type1"],"_"))[1:2],collapse="_")
  
  
  allele1 <- c(patientID,  sname,  hla,  x["HLA_A_type1"],  x["Allele_Freq_expected_HLA_type1"], x["Allele_Freq_Observed_HLA_type1"], x["HLALogR_Tumor_vs_Normal"], 
    x["HLA_type1copyNum_withBAF"], x["HLA_type1copyNum_withBAF_lower"], x["HLA_type1copyNum_withBAF_upper"], 
    x["HLA_type1copyNum_withBAFBin"], x["HLA_type1copyNum_withBAFBin_lower"], x["HLA_type1copyNum_withBAFBin_upper"], 
    x["HLA_type1copyNum_withoutBAF"], x["HLA_type1copyNum_withoutBAF_lower"], x["HLA_type1copyNum_withoutBAF_upper"], 
    x["HLA_type1copyNum_withoutBAFBin"], x["HLA_type1copyNum_withoutBAFBin_lower"], x["HLA_type1copyNum_withoutBAFBin_upper"], 
    x["PVal"], x["PVal_unique"], 
    x["minCoverageFilter_gene"], x["Total_Diff_Pos_between_alleles"], x["propSupportiveSites"], 
    allele_1_mismatch, allele_1_gaps,  x["Pass_HLA_Allele1_Diff_pos_Normal"])
  
  allele2 <- c(patientID,  sname,  hla,  x["HLA_A_type2"],  1-as.numeric(x["Allele_Freq_expected_HLA_type1"]), 
               1-as.numeric(x["Allele_Freq_Observed_HLA_type1"]), x["HLALogR_Tumor_vs_Normal"], 
               x["HLA_type2copyNum_withBAF"], x["HLA_type2copyNum_withBAF_lower"], x["HLA_type2copyNum_withBAF_upper"], 
               x["HLA_type2copyNum_withBAFBin"], x["HLA_type2copyNum_withBAFBin_lower"], x["HLA_type2copyNum_withBAFBin_upper"], 
               x["HLA_type2copyNum_withoutBAF"], x["HLA_type2copyNum_withoutBAF_lower"], x["HLA_type2copyNum_withoutBAF_upper"], 
               x["HLA_type2copyNum_withoutBAFBin"], x["HLA_type2copyNum_withoutBAFBin_lower"], x["HLA_type2copyNum_withoutBAFBin_upper"], 
               x["PVal"], x["PVal_unique"], 
               x["minCoverageFilter_gene"], x["Total_Diff_Pos_between_alleles"], x["propSupportiveSites"], 
               allele_1_mismatch, allele_1_gaps,  x["Pass_HLA_Allele1_Diff_pos_Normal"])
  
  return(list(allele1,allele2))
  
  }),recursive = FALSE))
  colnames(all)[1:3] <- c("mrn","sample","class")
  final_df <- data.frame(rbind(df,all)) %>% dplyr::distinct()
} else {
final_df <- df
}

write.table(final_df,file=paste(rdatadir,sname,"_hla_cn.txt", sep=''), sep="\t", row.names=FALSE, quote= FALSE)
