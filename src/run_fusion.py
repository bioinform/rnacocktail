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

def run_fusioncatcher(data_dir="", input="",  start=0, 
                  fusioncatcher=FUSIONCATCHER, fusioncatcher_opts="",  
                  sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT):


    logger.info("Running RNA fusion detection (FusionCatcher) for %s"%sample)
    if not os.path.exists(data_dir):
        logger.error("Aborting!")
        raise Exception("No data directory %s"%data_dir)


    work_fusioncatcher=os.path.join(workdir,"fusioncatcher",sample)
    create_dirs([work_fusioncatcher])
    fusioncatcher_log = os.path.join(work_fusioncatcher, "fusioncatcher.log")
    fusioncatcher_log_fd = open(fusioncatcher_log, "w")
    
    if nthreads>1:
        if "-p " not in fusioncatcher_opts:
            fusioncatcher_opts += " -p %d"%nthreads 
    msg = "Run FusionCatcher for %s"%sample
    command="%s -d %s -i %s --start %d -o %s" % (
        fusioncatcher, data_dir, input, start, work_fusioncatcher)
    command="bash -c \"%s\""%command        
    cmd = TimedExternalCmd(command, logger, raise_exception=True)
    retcode = cmd.run(cmd_log_fd_out=fusioncatcher_log_fd, cmd_log=fusioncatcher_log_fd, msg=msg, timeout=timeout)

    out_fusioncatcher=os.path.join(outdir,"fusioncatcher",sample)
    create_dirs([out_fusioncatcher])
    msg="Copy predictions to output directory for %s."%sample
    if os.path.exists("%s/final-list_candidate-fusion-genes.txt"%work_fusioncatcher):
        command = "cp %s/final-list_candidate-fusion-genes.txt %s/final-list_candidate-fusion-genes.txt"%(
                   work_fusioncatcher, out_fusioncatcher)
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=fusioncatcher_log_fd, cmd_log=fusioncatcher_log, msg=msg, timeout=timeout)   

    fusions = ""
    if os.path.exists("%s/final-list_candidate-fusion-genes.txt"%out_fusioncatcher):
        logger.info("FusionCatcher was successfull!")
        logger.info("Output fusions: %s/final-list_candidate-fusion-genes.txt"%out_fusioncatcher)
        fusions = "%s/final-list_candidate-fusion-genes.txt"%out_fusioncatcher
    else:            
        logger.info("FusionCatcher failed!")
    return fusions
    

def run_fusion(fusion_caller="FusionCatcher",
                  data_dir="", input="", start=0, 
                  fusioncatcher=FUSIONCATCHER, fusioncatcher_opts="",  
                  sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT, ignore_exceptions=False):
    fusions=""
    if fusion_caller.upper()=="FUSIONCATCHER":
        try:
            fusions=run_fusioncatcher(data_dir=data_dir, input=input, start=start, 
                      fusioncatcher=fusioncatcher, fusioncatcher_opts=fusioncatcher_opts,  
                      sample= sample, nthreads=nthreads, 
                      workdir=workdir, outdir=outdir, timeout=timeout)
        except Exception as excp:
            logger.info("FusionCatcher failed!")
            logger.error(excp)
            if not ignore_exceptions:
                raise Exception(excp)
    return fusions
    
    

