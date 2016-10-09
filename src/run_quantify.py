import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def run_salmon_smem(quantifier_idx=None,
                  seq_1="", seq_2="", seq_u="",
                  salmon_k=SALMON_SMEM_k, libtype="",
                  salmon_smem_opts="", salmon=SALMON,
                  start=0, sample= "", nthreads=1, unzip=False,
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running quantification (Salmon-SMEM) for %s"%sample)
    if not os.path.exists(quantifier_idx):
        logger.error("Aborting!")
        raise Exception("No Salmon FMD index directory %s"%quantifier_idx)
        
    if seq_1 and seq_2:
        for s1 in seq_1.split(","):
            if not os.path.exists(s1):
                logger.error("Aborting!")
                raise Exception("No Mate 1 sequence file %s"%s1)
                               
        for s2 in seq_2.split(","):
            if not os.path.exists(s2):
                logger.error("Aborting!")
                raise Exception("No Mate 2 sequence file %s"%s2)

        if unzip:
            seq_argument="-1 <(gunzip -c %s) -2 <(gunzip -c %s)"%(" ".join(seq_1.split(","))," ".join(seq_2.split(",")))
        else:
            if "," in seq_1:
                seq_1="<(cat %s)"%(" ".join(seq_1.split(",")))
            if "," in seq_2:
                seq_2="<(cat %s)"%(" ".join(seq_2.split(",")))
            seq_argument="-1 %s -2 %s"%(seq_1,seq_2)
    elif seq_u:
        if unzip:
            seq_argument="-r <(gunzip -c %s)"%(" ".join(seq_u.split(",")))
        elif "," in seq_u:
               seq_argument="-r <(cat %s)"%(" ".join(seq_u1.split(",")))
        else:
               seq_argument="-r %s"%(seq_u)
        for su in seq_u.split(","):
            if not os.path.exists(su):
                logger.error("Aborting!")
                raise Exception("No unpaired sequence file %s"%su)

    
    work_salmon_smem=os.path.join(workdir,"salmon_smem",sample)
    create_dirs([work_salmon_smem])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase Salmon-SMEM work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_salmon_smem)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1


    salmon_smem_log = os.path.join(work_salmon_smem, "salmon_smem.log")
    salmon_smem_log_fd = open(salmon_smem_log, "w")

    if "-p " not in salmon_smem_opts:
        salmon_smem_opts += " -p %d"%nthreads 

    salmon_smem_opts += " -k %d"%salmon_k 
    salmon_smem_opts += " -l %s"%libtype 

    msg = "Salmon-SMEM for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s quant -i %s %s %s -o %s" % (
            salmon, quantifier_idx, salmon_smem_opts, seq_argument,work_salmon_smem )
        command="bash -c \"%s\""%command
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=salmon_smem_log_fd, cmd_log=salmon_smem_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    

    out_salmon_smem=os.path.join(outdir,"salmon_smem",sample)
    create_dirs([out_salmon_smem])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/quant.sf"%work_salmon_smem):
            command = "cp %s/quant.sf %s/quant.sf"%(
                       work_salmon_smem, out_salmon_smem)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=salmon_smem_log_fd, cmd_log=salmon_smem_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    quant = ""
    if os.path.exists("%s/quant.sf"%out_salmon_smem):
        logger.info("Salmon-SMEM was successfull!")
        logger.info("Output expressions: %s/quant.sf"%out_salmon_smem)
        quant = "%s/quant.sf"%out_salmon_smem
    else:            
        logger.info("Salmon-SMEM was not successfull!")
    return quant

def run_quantify(quantifier="Salmon-SMEM", quantifier_idx=None,
                  seq_1="", seq_2="", seq_u="",
                  salmon_k=SALMON_SMEM_k, libtype="",
                  salmon_smem_opts="", salmon=SALMON,
                  start=0, sample= "", nthreads=1, unzip=False,
                  workdir=None, outdir=None, timeout=TIMEOUT):
    quant=""
    if quantifier.upper()=="SALMON-SMEM":
        quant=run_salmon_smem(quantifier_idx=quantifier_idx,
                      seq_1=seq_1, seq_2=seq_2, seq_u=seq_u,
                      salmon_k=salmon_k, libtype=libtype,
                      salmon_smem_opts=salmon_smem_opts, salmon=salmon,
                      start=start, sample= sample, nthreads=nthreads, unzip=unzip,
                      workdir=workdir, outdir=outdir, timeout=timeout)
    return quant