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

def run_lordec(kmer=23,
                  solid=3, long="", short="",
                  lordec=LORDEC, lordec_opts="",
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running long read error correction (LoRDEC) for %s"%sample)
    if not os.path.exists(long):
        logger.error("Aborting!")
        raise Exception("No long read sequence file %s"%long)

    if not os.path.exists(short):
        logger.error("Aborting!")
        raise Exception("No short read sequence file %s"%short)

    work_lordec=os.path.join(workdir,"lordec",sample)
    create_dirs([work_lordec])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase LoRDEC work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_lordec)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1

    lordec_log = os.path.join(work_lordec, "lordec.log")
    lordec_log_fd = open(lordec_log, "w")
    ksps = ""

    if "-T " not in lordec_opts:
        lordec_opts += " -T %d"%nthreads 

    msg = "LoRDEC for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s  -k %d -s %d -i %s -2 %s -O %s/../../../tmp -o %s/long_corrected.fa" % (
            lordec, lordec_opts, kmer, solid, long, short, work_lordec, work_lordec)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=lordec_log_fd, cmd_log=lordec_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    
    out_lordec=os.path.join(outdir,"lordec",sample)
    create_dirs([out_lordec])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/long_corrected.fa"%work_lordec):
            command = "cp %s/long_corrected.fa %s/long_corrected.fa"%(
                       work_lordec, out_lordec)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=lordec_log_fd, cmd_log=lordec_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    corrected = ""
    if os.path.exists("%s/long_corrected.fa"%out_lordec):
        logger.info("LoRDEC was successfull!")
        logger.info("Output corrected reads: %s/long_corrected.fa"%out_lordec)
        corrected = "%s/long_corrected.fa"%out_lordec
    else:            
        logger.info("LoRDEC failed!")
    return corrected

def run_lr_correct(long_corrector="LoRDEC", kmer=23,
                  solid=3, long="", short="",
                  lordec=LORDEC, lordec_opts="",
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT, ignore_exceptions=False):
    corrected=""
    if long_corrector.upper()=="LORDEC":
        try:
            corrected=run_lordec(kmer=kmer, solid=solid, long=long, short=short,
                      lordec=lordec, lordec_opts=lordec_opts,
                      start=start, sample= sample, nthreads=nthreads,
                      workdir=workdir, outdir=outdir, timeout=timeout)
        except Exception as excp:
            logger.info("LoRDEC failed!")
            logger.error(excp)
            if not ignore_exceptions:
                raise Exception(excp)
    return corrected