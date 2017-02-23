import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

def run_stringtie(alignment_bam="",ref_gtf="", 
                  stringtie_opts="", stringtie=STRINGTIE,
                  start=0, sample= "", nthreads=1,
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running transcriptome reconstruction (StringTie) for %s"%sample)
    if not os.path.exists(alignment_bam):
        logger.error("Aborting!")
        raise Exception("No input alignment BAM file %s"%alignment_bam)
        
    work_stringtie="%s/stringtie/%s/"%(workdir,sample)
    create_dirs([work_stringtie])
    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase StringTie work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_stringtie)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1
    stringtie_log = os.path.join(work_stringtie, "stringtie.log")
    stringtie_log_fd = open(stringtie_log, "w")

    if ref_gtf:
        if not os.path.exists(ref_gtf):
            logger.error("Aborting!")
            raise Exception("No reference GTF file %s"%ref_gtf)

    if ref_gtf:
        stringtie_opts += " -G %s"%ref_gtf
    if "-p " not in stringtie_opts:
        stringtie_opts += " -p %d"%nthreads 

    msg = "StringTie for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s %s -o %s/transcripts.gtf -A %s/gene_abund.tab -v" % (
            stringtie, alignment_bam, stringtie_opts, work_stringtie, work_stringtie)
        command="bash -c \"%s\""%command
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=stringtie_log_fd, cmd_log=stringtie_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    
    out_stringtie=os.path.join(outdir,"stringtie",sample)
    create_dirs([out_stringtie])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/transcripts.gtf"%work_stringtie) and \
           os.path.exists("%s/gene_abund.tab"%work_stringtie):
            command = "cp %s/transcripts.gtf %s/transcripts.gtf"%(
                       work_stringtie, out_stringtie)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=stringtie_log_fd, cmd_log=stringtie_log, msg=msg, timeout=timeout)   
            
            command = "cp %s/gene_abund.tab %s/gene_abund.tab"%(
                       work_stringtie, out_stringtie)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=stringtie_log_fd, cmd_log=stringtie_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    transcripts = ""
    abundances = ""
    if os.path.exists("%s/transcripts.gtf"%out_stringtie) and \
       os.path.exists("%s/gene_abund.tab"%out_stringtie):
        logger.info("StringTie was successfull!")
        logger.info("Output isoforms: %s/transcripts.gtf"%out_stringtie)
        logger.info("Output expressions: %s/gene_abund.tab"%out_stringtie)
        transcripts = "%s/transcripts.gtf"%out_stringtie   
        abundances = "%s/gene_abund.tab"%out_stringtie   
    else:            
        logger.info("StringTie failed!")
    return transcripts,abundances

def run_reconstruct(reconstructor="StringTie", alignment_bam="",
                  ref_gtf="", 
                  stringtie_opts="", stringtie=STRINGTIE,
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT, ignore_exceptions=False):
    transcripts = ""
    abundances = ""
    if reconstructor.upper()=="STRINGTIE":
        try:
            transcripts,abundances=run_stringtie(alignment_bam=alignment_bam,
                          ref_gtf=ref_gtf, 
                          stringtie_opts=stringtie_opts, stringtie=stringtie,
                          start=start, sample= sample, nthreads=nthreads,
                          workdir=workdir, outdir=outdir, timeout=timeout)
        except Exception as excp:
            logger.info("StringTie failed!")
            if not ignore_exceptions:
                raise Exception(excp)
            else:
                logger.error(excp)
    return transcripts,abundances