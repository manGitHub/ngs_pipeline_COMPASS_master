import sys

# for strelka indel and snp vcfs
# adds AD to FORMAT and sample columns
# adds VAF_1 and VAF_2 to INFO columns

# usage:
# python3 <sample>.strelka.indels.raw.vcf <sample>.strelka.indels.raw.AD.vcf

fileN = sys.argv[1]

""" Check for snp or indel in strelka file """
def check_if_string_in_file(file_name, string_to_search):
    # Open the file in read only mode
    with open(file_name, 'r') as read_obj:
        # Read all lines in the file one by one
        for line in read_obj:
            # For each line, check if line contains the string
            if string_to_search in line:
                return True
    return False

""" snp strelka file """
def vaf_snp(file_name):
    exitFlag = False
    with open(file_name) as myfile:
        info = []
        notInfo = []
        chrom = []
        dat = []
        format = []
        for each in myfile:                                                     # for each line in file
            if each.startswith('##'):                                           # subset commented section
                if each.startswith('##INFO'):                                   # subset info lines, will add vaf later
                    info.append(each.strip())
                elif each.startswith('##FORMAT'):
                    format.append(each.strip())
                else:                                                           # if not info lines, subset
                    notInfo.append(each.strip())
            if each.startswith('#CHROM'):                                       # pull out table header
                chrom.append(each.strip())

            if not each.startswith('#'):                                        # subset table and parse ACGTs, use tier 1

                # error handling correct vcf formatting
                if each.split('\t')[8]!='DP:FDP:SDP:SUBDP:AU:CU:GU:TU':
                    exitFlag = True

                ## First Sample
                A_1 = float(each.split('\t')[9].split(':')[4].split(',')[0])
                C_1 = float(each.split('\t')[9].split(':')[5].split(',')[0])
                G_1 = float(each.split('\t')[9].split(':')[6].split(',')[0])
                T_1 = float(each.split('\t')[9].split(':')[7].split(',')[0])

                ## Second Sample
                As = float(each.split('\t')[10].split(':')[4].split(',')[0])
                Cs = float(each.split('\t')[10].split(':')[5].split(',')[0])
                Gs = float(each.split('\t')[10].split(':')[6].split(',')[0])
                Ts = float(each.split('\t')[10].split(':')[7].split(',')[0])

                # First sample: DEPTH_1
                # Second sample: DEPTH
                REF = each.split('\t')[3]                                      # Find REF variant, set appropriate value for depth
                if REF == 'A':
                    DEPTH = As
                    DEPTH_1 = A_1
                elif REF == 'C':
                    DEPTH = Cs
                    DEPTH_1 = C_1
                elif REF == 'G':
                    DEPTH = Gs
                    DEPTH_1 = G_1
                elif REF == 'T':
                    DEPTH = Ts
                    DEPTH_1 = T_1

                # First sample: VAF_1, AD_1
                # Second sample: VAF, AD
                ALT = each.split('\t')[4]                                       # find ALT variant, calc VAF
                if ALT == 'A':
                    VAF_1 = round(A_1 / (A_1 + DEPTH_1), 2) if (A_1 + DEPTH_1)>0 else 0
                    VAF = round(As / (As + DEPTH), 2) if (As + DEPTH)>0 else 0
                    AD_1 = A_1
                    AD = As
                elif ALT == 'C':
                    VAF_1 = round(C_1 / (C_1 + DEPTH_1), 2) if (C_1 + DEPTH_1)>0 else 0
                    VAF = round(Cs / (Cs + DEPTH), 2) if (Cs + DEPTH)>0 else 0
                    AD_1 = C_1
                    AD = Cs
                elif ALT == 'G':
                    VAF_1 = round(G_1 / (G_1 + DEPTH_1), 2) if (G_1 + DEPTH_1)>0 else 0
                    VAF = round(Gs / (Gs + DEPTH), 2) if (Gs + DEPTH)>0 else 0
                    AD_1 = G_1
                    AD = Gs
                elif ALT == 'T':
                    VAF_1 = round(T_1 / (T_1 + DEPTH_1), 2) if (T_1 + DEPTH_1)>0 else 0
                    VAF = round(Ts / (Ts + DEPTH), 2) if (Ts + DEPTH)>0 else 0
                    AD_1 = T_1
                    AD = Ts

                # add VAF to INFO column
                sepEach = each.split('\t')
                sepEach[7] = sepEach[7] + ';VAF_1=' + str(VAF_1) + ';VAF_2=' + str(VAF)                                 # new info column with vaf

                # add AD to FORMAT and sample columns
                sepEach[8] = sepEach[8] + ':AD'
                sepEach[9] = sepEach[9] + ':' + str(int(DEPTH_1)) + ',' + str(int(AD_1))                                # first Sample
                sepEach[10] = sepEach[10].strip() + ':' + str(int(DEPTH)) + ',' + str(int(AD))                          # second Sample

                eachStr = '\t'
                each = eachStr.join(sepEach)                                                                            # redefine orig each variable
                dat.append(each.strip())

        info.append('##INFO=<ID=VAF_1,Number=1,Type=Float,Description="Variant Allele Frequency for First Sample">')    # add info tag vaf to header
        info.append('##INFO=<ID=VAF_2,Number=1,Type=Float,Description="Variant Allele Frequency for Second Sample">')
        format.append('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">')
        together = notInfo + format + info + chrom + dat                                                                         # put vcf file back together
        if exitFlag:
            sys.exit('Error: vcf FORMAT field is not correct, should be DP:FDP:SDP:SUBDP:AU:CU:GU:TU')
        else:
            outFile = open(sys.argv[2],'w')
            for x in together:
                outFile.write(x.strip() + '\n')
        outFile.close()

""" indel strelka file """
def vaf_indel(file_name):
    exitFlag = False
    with open(file_name) as myfile:
        info = []
        notInfo = []
        chrom = []
        dat = []
        format = []
        for each in myfile:                                                     # for each line in file
            if each.startswith('##'):                                           # subset commented section
                if each.startswith('##INFO'):                                   # subset info and format lines, will add vaf and AD later
                    info.append(each.strip())
                elif each.startswith('##FORMAT'):
                    format.append(each.strip())
                else:                                                           # if not info or format header, subset
                    notInfo.append(each.strip())
            if each.startswith('#CHROM'):                                       # pull out table header
                chrom.append(each.strip())
            if not each.startswith('#'):                                        # subset table and parse TAR and TIR, use tier 1

                # error handling correct vcf formatting
                if each.split('\t')[8]!='DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50':
                    exitFlag = True

                # First sample
                REF_1 = float(each.split('\t')[9].split(':')[2].split(',')[0])   # REF counts (TAR)
                ALT_1 = float(each.split('\t')[9].split(':')[3].split(',')[0])   # ALT counts (TIR)
                VAF_1 = round(ALT_1 / (ALT_1 + REF_1), 2) if (ALT_1 + REF_1)>0 else 0

                # Second sample
                REF = float(each.split('\t')[10].split(':')[2].split(',')[0])   # REF counts (TAR)
                ALT = float(each.split('\t')[10].split(':')[3].split(',')[0])   # ALT counts (TIR)
                VAF = round(ALT / (ALT + REF), 2) if (ALT + REF)>0 else 0

                # add vaf to INFO column
                sepEach = each.split('\t')
                sepEach[7] = sepEach[7] + ';VAF_1=' + str(VAF_1) + ';VAF_2=' + str(VAF)                    # new info column with vaf

                # add AD to FORMAT column and sample columns
                sepEach[8] = sepEach[8] + ':AD'                                 # FORMAT
                sepEach[9] = sepEach[9] + ':' + sepEach[9].split(':')[2].split(',')[0] + ',' + sepEach[9].split(':')[3].split(',')[0]       # first sample
                sepEach[10] = sepEach[10].strip() + ':' + sepEach[10].split(':')[2].split(',')[0] + ',' + sepEach[10].split(':')[3].split(',')[0]   # second sample

                eachStr = '\t'
                each = eachStr.join(sepEach)                                    # redefine orig each variable
                dat.append(each.strip())

        info.append('##INFO=<ID=VAF_1,Number=1,Type=Float,Description="Variant Allele Frequency for First Sample">')   # add info tag vaf to header
        info.append('##INFO=<ID=VAF_2,Number=1,Type=Float,Description="Variant Allele Frequency for Second Sample">')   # add info tag vaf to header
        format.append('##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">')
        together = notInfo + format + info + chrom + dat                                                     # put vcf file back together
        if exitFlag:
            sys.exit('Error: vcf FORMAT field is not correct, should be DP:DP2:TAR:TIR:TOR:DP50:FDP50:SUBDP50:BCN50')
        else:
            outFile = open(sys.argv[2],'w')
            for x in together:
                outFile.write(x.strip() + '\n')
        outFile.close()

""" Check type of file, run appropriate function """
if check_if_string_in_file(fileN,'##content=strelka somatic snv calls'):
    vaf_snp(fileN)
elif check_if_string_in_file(fileN,'##content=strelka somatic indel calls'):
    vaf_indel(fileN)
