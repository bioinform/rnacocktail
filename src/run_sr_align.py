import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def run_hisat2(align_idx=None,
                  seq_1="", seq_2="", seq_u="",
                  seq_sra="", ref_gtf="", 
                  hisat2_opts="", hisat2=HISAT2, hisat2_sps=HISAT2_SPS,
                  samtools=SAMTOOLS,
                  start=0, sample= "", nthreads=1,
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running alignment (HISAT2) for %s"%sample)
    if not os.path.exists(align_idx+".1.ht2"):
        logger.error("Aborting!")
        raise Exception("No HISAT index file %s.1.ht2"%align_idx)
        
    if seq_1 and seq_2:
        for s1 in seq_1.split(","):
            if not os.path.exists(s1):
                logger.error("Aborting!")
                raise Exception("No Mate 1 sequence file %s"%s1)
        for s2 in seq_2.split(","):
            if not os.path.exists(s2):
                logger.error("Aborting!")
                raise Exception("No Mate 2 sequence file %s"%s2)
        seq_argument="-1 %s -2 %s"%(seq_1,seq_2)
    elif seq_u:
        seq_argument="-U %s"%(seq_u)
        for su in seq_u.split(","):
            if not os.path.exists(su):
                logger.error("Aborting!")
                raise Exception("No unpaired sequence file %s"%su)

    elif seq_sra:
        seq_argument="--sra-acc %s"%(seq_sra)
        for sr in seq_sra.split(","):
            if not os.path.exists(sr):
                logger.error("Aborting!")
                raise Exception("No sra sequence file %s"%sr)


    work_hisat2=os.path.join(workdir,"hisat2",sample)
    create_dirs([work_hisat2])
 
    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase HISAT2 work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_hisat2)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1

    hisat2_log = os.path.join(work_hisat2, "hisat2.log")
    hisat2_log_fd = open(hisat2_log, "w")
    
    ksps = ""
    msg = "Prepare known-splicesites for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if ref_gtf:
            if not os.path.exists(ref_gtf):
                logger.error("Aborting!")
                raise Exception("No reference GTF file %s"%ref_gtf)
            else:
                ksps =  ref_gtf.strip() + "known-splicesite.txt"
                if os.path.exists(ksps):
                    logger.info("Will use the precomputed %s as --known-splicesite-infile for HISAT2"%ksps)
                else:
                    msg="compute --known-splicesite-infile for HISAT2"
                    ksps =  os.path.join(work_hisat2, "known-splicesite.txt")
                    ksps_fd = open(ksps, "w")
                
                    command="%s %s" % (hisat2_sps,ref_gtf)
                    command="bash -c \"%s\""%command
                    cmd = TimedExternalCmd(command, logger, raise_exception=True)
                    retcode = cmd.run(cmd_log_fd_out=ksps_fd, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    

    
    if "--dta " not in hisat2_opts:
        hisat2_opts += " --dta"
    if "--rg-id " not in hisat2_opts:
        hisat2_opts += " --rg-id hisat2"
    if "--rg " not in hisat2_opts:
        hisat2_opts += " --rg SM:%s"%sample
    if "--threads " not in hisat2_opts:
        hisat2_opts += " --threads %d"%nthreads 
    if ksps:
        hisat2_opts += " --known-splicesite-infile %s"%ksps

    msg = "HISAT2 for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s  -x %s %s -S %s/alignments.sam --novel-splicesite-outfile %s/splicesites.tab" % (
            hisat2, hisat2_opts, align_idx, seq_argument,work_hisat2, work_hisat2 )
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=hisat2_log_fd, cmd_log=hisat2_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "converting SAM to BAM for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s view -Su %s/alignments.sam -@ %d -o %s/alignments.bam" % (
            samtools, work_hisat2, nthreads, work_hisat2)
        command="bash -c \"%s\""%command       
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=hisat2_log_fd, cmd_log=hisat2_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "sorting BAM for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s sort  -@ %d %s/alignments.bam %s/alignments.sorted " % (
            samtools, nthreads, work_hisat2, work_hisat2)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=hisat2_log_fd, cmd_log=hisat2_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    


    msg = "Converting junctions to BED for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="hisat2_jun2bed.py %s/splicesites.tab %s/splicesites.bed " % (
             work_hisat2, work_hisat2)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=hisat2_log_fd, cmd_log=hisat2_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

#     msg = "Clean temp alignment files for %s"%sample
#     if start<=step:
#         logger.info("--------------------------STEP %s--------------------------"%step)
#         command="rm %s/alignments.sam %s/alignments.bam" % (work_hisat2, work_hisat2)
#         command="bash -c \"%s\""%command    
#         cmd = TimedExternalCmd(command, logger, raise_exception=True)
#         retcode = cmd.run(cmd_log_fd_out=hisat2_log_fd, cmd_log=hisat2_log, msg=msg, timeout=timeout)
#     else:
#         logger.info("Skipping step %d: %s"%(step,msg))
#     step+=1


    out_hisat2=os.path.join(outdir,"hisat2",sample)
    create_dirs([out_hisat2])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/alignments.sorted.bam"%work_hisat2):
            command = "cp %s/alignments.sorted.bam %s/alignments.sorted.bam"%(
                       work_hisat2, out_hisat2)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=hisat2_log_fd, cmd_log=hisat2_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1



    alignments_bam = ""
    junctions_tab = ""
    junctions_bed = ""
    if os.path.exists("%s/alignments.sorted.bam"%out_hisat2):
        logger.info("HISAT2 was successfull!")
        logger.info("Output alignment: %s/alignments.sorted.bam"%out_hisat2)
        logger.info("Output junction tab: %s/splicesites.tab"%out_hisat2)
        logger.info("Output junction bed: %s/splicesites.bed"%out_hisat2)
        alignments_bam = "%s/alignments.sorted.bam"%out_hisat2   
        junctions_tab = "%s/splicesites.tab"%out_hisat2   
        junctions_bed = "%s/splicesites.bed"%out_hisat2   
    else:            
        logger.info("HISAT2 was not successfull!")
    return alignments_bam,junctions_tab,junctions_bed

def run_sr_align(sr_aligner="HISAT2", align_idx=None,
                  seq_1="", seq_2="", seq_u="",
                  seq_sra="", ref_gtf="", 
                  hisat2_opts="", hisat2=HISAT2, hisat2_sps=HISAT2_SPS,
                  samtools=SAMTOOLS,
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT):
    alignments_bam=""
    junctions_tab = ""
    junctions_bed = ""
    if sr_aligner.upper()=="HISAT2":
        alignments_bam=run_hisat2(align_idx=align_idx,
                      seq_1=seq_1, seq_2=seq_2, seq_u=seq_u,
                      seq_sra=seq_sra, ref_gtf=ref_gtf, 
                      hisat2_opts=hisat2_opts, hisat2=hisat2, hisat2_sps=hisat2_sps,
                      samtools=samtools,
                      start=start, sample= sample, nthreads=nthreads,
                      workdir=workdir, outdir=outdir, timeout=timeout)
    return alignments_bam,junctions_tab,junctions_bed