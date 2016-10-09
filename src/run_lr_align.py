import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def run_starlong(long="", 
                  genome_dir="", ref_gtf="",
                  starlong=STARLONG, sam2psl=SAM2PSL,samtools=SAMTOOLS,
                  starlong_opts="",
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running long read alignment (STARlong) for %s"%sample)
    if not os.path.exists(genome_dir+"SAindex"):
        logger.error("Aborting!")
        raise Exception("No SAindex directory in %s"%genome_dir)
        
    if long:
        if not os.path.exists(long):
            logger.error("Aborting!")
            raise Exception("No long read sequence file %s"%long)

    work_starlong=os.path.join(workdir,"starlong",sample)
    create_dirs([work_starlong])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase STARlong work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_starlong)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1
    
    starlong_log = os.path.join(work_starlong, "starlong.log")
    starlong_log_fd = open(starlong_log, "w")
    
    
    
    if ref_gtf:
        if not os.path.exists(ref_gtf):
            logger.error("Aborting!")
            raise Exception("No reference GTF file %s"%ref_gtf)    

    if "--outSAMattrRGline" not in starlong_opts:
        starlong_opts += " --outSAMattrRGline ID:STARlong SM:%s"%sample
    if "--runThreadN " not in starlong_opts:
        starlong_opts += " --runThreadN %d"%nthreads 
    if ref_gtf:
        starlong_opts += " --sjdbGTFfile %s"%ref_gtf 
    for k,v in STARLONG_DEFAULTS.iteritems():
        if k not in starlong_opts:
            starlong_opts += " --%s %s"%(k,v) 


    msg = "STARlong for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s --runMode alignReads %s --genomeDir %s  --readFilesIn %s  --outFileNamePrefix %s/" % (
            starlong, starlong_opts, genome_dir, long, work_starlong )
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=starlong_log_fd, cmd_log=starlong_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    
    msg = "converting SAM to PSL for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s -i %s/Aligned.out.sam -o %s/Aligned.out.psl" % (
            sam2psl, work_starlong, work_starlong)
        command="bash -c \"%s\""%command       
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=starlong_log_fd, cmd_log=starlong_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "converting SAM to BAM for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s view -Su %s/Aligned.out.sam -@ %d -o %s/Aligned.out.bam" % (
            samtools, work_starlong, nthreads, work_starlong)
        command="bash -c \"%s\""%command       
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=starlong_log_fd, cmd_log=starlong_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

# 
#     msg = "Clean temp alignment files for %s"%sample
#     if start<=step:
#         logger.info("--------------------------STEP %s--------------------------"%step)
#         command="rm %s/Aligned.out.sam" % (work_starlong)
#         command="bash -c \"%s\""%command    
#         cmd = TimedExternalCmd(command, logger, raise_exception=True)
#         retcode = cmd.run(cmd_log_fd_out=starlong_log_fd, cmd_log=starlong_log, msg=msg, timeout=timeout)
#     else:
#         logger.info("Skipping step %d: %s"%(step,msg))
#     step+=1


    out_starlong=os.path.join(outdir,"starlong",sample)
    create_dirs([out_starlong])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/Aligned.out.psl"%work_starlong):
            command = "cp %s/Aligned.out.psl %s/Aligned.out.psl"%(
                       work_starlong, out_starlong)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=starlong_log_fd, cmd_log=starlong_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    alignments_psl = ""
    if os.path.exists("%s/Aligned.out.psl"%out_starlong):
        logger.info("STARlong was successfull!")
        logger.info("Output alignment: %s/Aligned.out.psl"%out_starlong)
        alignments_psl = "%s/Aligned.out.psl"%out_starlong
    else:            
        logger.info("STARlong was not successfull!")
    return alignments_psl

def run_lr_align(long_aligner="STARlong", long="", 
                  genome_dir="", ref_gtf="",
                  starlong=STARLONG, sam2psl=SAM2PSL, samtools=SAMTOOLS,
                  starlong_opts="",
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT):
    alignments_psl=""
    if long_aligner.upper()=="STARLONG":
        alignments_bam=run_starlong(genome_dir=genome_dir, ref_gtf=ref_gtf,
                      long=long, starlong=starlong, sam2psl=sam2psl, samtools=samtools,
                      starlong_opts=starlong_opts,
                      start=start, sample= sample, nthreads=nthreads,
                      workdir=workdir, outdir=outdir, timeout=timeout) 
    return alignments_psl