import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def run_oases(assmebly_hash=DNV_HASH,
                seq_1="", seq_2="", seq_u="", seq_i="",
                file_format=DNV_FORMAT, read_type=DNV_READTYPE,
                oases=OASES, velvetg=VELVETG, velveth=VELVETH,
                oases_opts="", velvetg_opts="", velveth_opts="",
                start=0, sample= "", nthreads=1,
                workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running de novo assembly (OASES) for %s"%sample)

    if seq_1 and seq_2:
        for s1 in seq_1.split(","):
            if not os.path.exists(s1):
                logger.error("Aborting!")
                raise Exception("No Mate 1 sequence file %s"%s1)
        for s2 in seq_2.split(","):
            if not os.path.exists(s2):
                logger.error("Aborting!")
                raise Exception("No Mate 2 sequence file %s"%s2)
        seq_argument="-separate %s %s"%(seq_1,seq_2)
    elif seq_u:
        seq_argument=seq_u
        for su in seq_u.split(","):
            if not os.path.exists(su):
                logger.error("Aborting!")
                raise Exception("No unpaired sequence file %s"%su)

    elif seq_i:
        seq_argument=seq_i
        for sr in seq_i.split(","):
            if not os.path.exists(seq_i):
                logger.error("Aborting!")
                raise Exception("No sra sequence file %s"%sr)

    work_oases=os.path.join(workdir,"oases",sample)
    create_dirs([work_oases])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase Oases work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_oases)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg, timeout=timeout)
    step+=1

    oases_log = os.path.join(work_oases, "oases.log")
    oases_log_fd = open(oases_log, "w")

    
    seq_argument="-%s -%s %s "%(file_format,read_type,seq_argument)

    msg = "velveth for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s %d  %s %s" % (
            velveth, work_oases, assmebly_hash, velveth_opts, seq_argument)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True, env_dict={"OMP_NUM_THREADS":str(nthreads)})
        retcode = cmd.run(cmd_log_fd_out=oases_log_fd, cmd_log=oases_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    
    
    msg = "velvetg for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s %s -read_trkg yes " % (
            velvetg, work_oases, velvetg_opts)
        command="bash -c \"%s\""%command       
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=oases_log_fd, cmd_log=oases_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "oases for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s %s " % (
            oases, work_oases, oases_opts)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=oases_log_fd, cmd_log=oases_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    out_oases=os.path.join(outdir,"oases",sample)
    create_dirs([out_oases])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/transcripts.fa"%work_oases):
            command = "cp %s/transcripts.fa %s/transcripts.fa"%(
                       work_oases, out_oases)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=oases_log_fd, cmd_log=oases_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    
    transcripts = ""
    if os.path.exists("%s/transcripts.fa"%out_oases):
        logger.info("Oases was successfull!")
        logger.info("Output transcripts: %s/transcripts.fa"%out_oases)
        transcripts = "%s/transcripts.fa"%out_oases
    else:            
        logger.info("Oases was not successfull!")
    return transcripts

def run_dnv_assemebly(assembler="Oases", assmebly_hash=DNV_HASH,
                      seq_1="", seq_2="", seq_u="", seq_i="",
                      file_format=DNV_FORMAT, read_type=DNV_READTYPE, 
                      oases=OASES, velvetg=VELVETG, velveth=VELVETH,
                      oases_opts="", velvetg_opts="", velveth_opts="",
                      start=0, sample= "", nthreads=1,
                      workdir=None, outdir=None, timeout=TIMEOUT):
    transcripts=""
    if assembler.upper()=="OASES":
        transcripts=run_oases(assmebly_hash=assmebly_hash,
                      seq_1=seq_1, seq_2=seq_2, seq_u=seq_u, seq_i=seq_i,
                      file_format=file_format, read_type=read_type, 
                      oases=oases, velvetg=velvetg, velveth=velveth,
                      oases_opts=oases_opts, velvetg_opts=velvetg_opts, velveth_opts=velveth_opts,
                      start=start, sample= sample, nthreads=nthreads,
                      workdir=workdir, outdir=outdir, timeout=timeout)
        
                      
    return transcripts