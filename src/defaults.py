MODES = set(["align", "reconstruct", "denovo",
             "quantify", "diff", "long_correct", "long_align",
             "long_reconstruct", "variant", "editing", "fusion"])
SR_ALIGNERS = set(["HISAT2"])
RECONSTRUCTORS = set(["StringTie"])
QUANTIFIERS = set(["Salmon-SMEM"])
DIFFS = set(["DESeq2"])
DNV_ASSEMBLERS = set(["Oases"])
LR_CORRECTORS = set(["LoRDEC"])
LR_ALIGNERS= set(["STARlong"])
LR_RECONSTRUCTORS= set(["IDP"])
variant_caller= set(["GATK"])
editing_caller= set(["GIRMI"])
fusion_caller= set(["FusionCatcher"])
TIMEOUT = 10000000  # in seconds


SALMON_SMEM_k = 19
DESeq2_MINCNT = 2
DESeq2_ALPHA = 0.05
DNV_HASH = 25
DNV_FORMAT = "fasta"
DNV_READTYPE = "short"
STARLONG_DEFAULTS = {"outSAMattributes": "NH HI NM MD", "readNameSeparator": "space",
                     "outFilterMultimapScoreRange": "1", "outFilterMismatchNmax": "2000",
                     "scoreGapNoncan": "-20", "scoreGapGCAG":"-4", "scoreGapATAC":"-8",
                     "scoreDelOpen": "-1", "scoreDelBase": "-1", "scoreInsOpen": "-1", "scoreInsBase": "-1", 
                     "alignEndsType": "Local", "seedSearchStartLmax": "50", "seedPerReadNmax": "100000", 
                     "seedPerWindowNmax": "1000", "alignTranscriptsPerReadNmax": "100000", 
                     "alignTranscriptsPerWindowNmax": "10000"}


GATK_SN_RF = "ReassignOneMappingQuality"
GATK_SN_RMQF = 255
GATK_SN_RMQT = 60
GATK_SN_OPT = (("-rf %s " % GATK_SN_RF) if GATK_SN_RF else "") + \
              (("-RMQF %s " % GATK_SN_RMQF) if GATK_SN_RMQF else "") + \
              (("-RMQT %s " % GATK_SN_RMQT) if GATK_SN_RMQT else "") + "-U ALLOW_N_CIGAR_READS"

GATK_HC_STANDCALLCONF = 20.0
GATK_HC_STANDEMITCONF = 20.0
GATK_HC_OPT = (("-stand_call_conf %f " % GATK_HC_STANDCALLCONF) if GATK_HC_STANDCALLCONF else "") + \
              (("-stand_emit_conf %f " % GATK_HC_STANDEMITCONF) if GATK_HC_STANDEMITCONF else "") + \
              "-dontUseSoftClippedBases"


GATK_VF_WINDOW = 35
GATK_VF_CLUSTER = 3
GATK_VF_FSMIN = 30.0
GATK_VF_QDMAX = 2.0
GATK_VF_OPT = (("-window %d " % GATK_VF_WINDOW) if GATK_VF_WINDOW else "") + \
              (("-cluster %d " % GATK_VF_CLUSTER) if GATK_VF_CLUSTER else "") + \
              (("-filterName FS -filter 'FS > %f' " % GATK_VF_FSMIN) if GATK_VF_FSMIN else "") + \
              (("-filterName QD -filter 'QD < %f' " % GATK_VF_QDMAX) if GATK_VF_QDMAX else "") 

JAVA_XMS = "-Xms1g"
JAVA_XMG = "-Xmx5g"
JAVA_OPT= "%s %s"%(JAVA_XMS,JAVA_XMG)


HISAT2 = "hisat2"
HISAT2_SPS = "hisat2_extract_splice_sites.py"
SAMTOOLS = "samtools"
STRINGTIE = "stringtie"
SALMON = "salmon"
R_CMD = "R"
FEATURECOUNTS = "featureCounts"
VELVETG = "velvetg"
VELVETH = "velveth"
OASES = "oases"
LORDEC = "lordec-correct"
STARLONG = "STARlong"
SAM2PSL = "sam2psl.py"
IDP = "runIDP.py"
PICARD = "picard.jar"
GATK = "GenomeAnalysisTK.jar"
JAVA = "java"
GIREMI = "giremi"
HTSLIB = ""
FUSIONCATCHER= "fusioncatcher"
