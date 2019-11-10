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

def run_gatk(alignment="", ref_genome="", knownsites="",
                  picard=PICARD, gatk=GATK,                  
                  java=JAVA, java_opts="",
                  CleanSam=False, no_BaseRecalibrator=False ,               
                  AddOrReplaceReadGroups_opts="",  MarkDuplicates_opts="",
                  SplitNCigarReads_opts="",
                  BaseRecalibrator_opts="",
                  ApplyBQSR_opts="",  HaplotypeCaller_opts="",
                  VariantFiltration_opts="",  
                  start=0, sample= "", nthreads=1,
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running variant calling (GATK) for %s"%sample)
    if not os.path.exists(alignment):
        logger.error("Aborting!")
        raise Exception("No alignment file %s"%alignment)
    if not os.path.exists(ref_genome):
        logger.error("Aborting!")
        raise Exception("No reference genome FASTA file %s"%ref_genome)


    work_gatk=os.path.join(workdir,"gatk",sample)
    create_dirs([work_gatk])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase GATK work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_gatk)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1

    gatk_log = os.path.join(work_gatk, "gatk.log")
    gatk_log_fd = open(gatk_log, "w")
    

    if "SO=" not in AddOrReplaceReadGroups_opts:
        AddOrReplaceReadGroups_opts += " SO=coordinate"
    if "RGLB=" not in AddOrReplaceReadGroups_opts:
        AddOrReplaceReadGroups_opts += " RGLB=lib1"
    if "RGPL=" not in AddOrReplaceReadGroups_opts:
        AddOrReplaceReadGroups_opts += " RGPL=illumina"
    if "RGPU=" not in AddOrReplaceReadGroups_opts:
        AddOrReplaceReadGroups_opts += " RGPU=unit1"
    if "RGSM=" not in AddOrReplaceReadGroups_opts:
        AddOrReplaceReadGroups_opts += " RGSM=%s"%sample

    if "CREATE_INDEX=" not in MarkDuplicates_opts:
        MarkDuplicates_opts += " CREATE_INDEX=true"
    if "VALIDATION_STRINGENCY=" not in MarkDuplicates_opts:
        MarkDuplicates_opts += " VALIDATION_STRINGENCY=SILENT"

    if knownsites:    
        if not os.path.exists(knownsites):
            logger.error("Aborting!")
            raise Exception("No VCF knownsites file %s"%knownsites)
        if "--known-sites " not in BaseRecalibrator_opts:
            BaseRecalibrator_opts += " --known-sites %s"%knownsites



    if "--dont-use-soft-clipped-bases " not in HaplotypeCaller_opts:
        HaplotypeCaller_opts += " --dont-use-soft-clipped-bases"
    if "-stand-call-conf " not in HaplotypeCaller_opts:
        HaplotypeCaller_opts += " -stand-call-conf %f"%GATK_HC_STANDCALLCONF

    if "-window " not in VariantFiltration_opts:
        VariantFiltration_opts += " -window %d"%GATK_VF_WINDOW
    if "-cluster " not in VariantFiltration_opts:
        VariantFiltration_opts += " -cluster %d"%GATK_VF_CLUSTER
    if "--filter-name FS " not in VariantFiltration_opts:
        VariantFiltration_opts += " --filter-name FS -filter 'FS > %f'"%GATK_VF_FSMIN
    if "--filter-name QD " not in VariantFiltration_opts:
        VariantFiltration_opts += " --filter-name QD -filter 'QD < %f'"%GATK_VF_QDMAX

    if "-Xms" not in java_opts:
        java_opts += " %s"%JAVA_XMS
    if "-Xmx" not in java_opts:
        java_opts += " %s"%JAVA_XMG
    if "-Djava.io.tmpdir" not in java_opts:
        java_opts += " -Djava.io.tmpdir=%s/javatmp/"%(work_gatk)
        create_dirs(["%s/javatmp/"%(work_gatk)])

    msg = "picard CleanSam for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if CleanSam:
            command="%s %s -cp %s picard.cmdline.PicardCommandLine CleanSam -I %s -O %s/alignments_clean.bam" % (
                java, java_opts, picard, alignment,work_gatk )
            command="bash -c \"%s\""%command      
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
            alignment="%s/alignments_clean.bam"%work_gatk
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    msg = "picard AddOrReplaceReadGroups for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s -cp %s picard.cmdline.PicardCommandLine AddOrReplaceReadGroups -I %s -O %s/rg_added_sorted.bam %s" % (
            java, java_opts, picard, alignment,work_gatk,AddOrReplaceReadGroups_opts)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    msg = "picard MarkDuplicates for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s -cp %s picard.cmdline.PicardCommandLine MarkDuplicates -I %s/rg_added_sorted.bam -O %s/dedupped.bam %s -M %s/output.metrics" % (
            java, java_opts, picard, work_gatk,work_gatk,MarkDuplicates_opts,work_gatk)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    msg = "GATK SplitNCigarReads for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s -jar %s SplitNCigarReads -R %s -I %s/dedupped.bam -O %s/split.bam %s" % (
            java, java_opts, gatk, ref_genome,work_gatk,work_gatk,SplitNCigarReads_opts)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    split_bam="%s/split.bam"%work_gatk

    if not no_BaseRecalibrator:
        msg = "GATK BaseRecalibrator for %s"%sample
        if start<=step:
            logger.info("--------------------------STEP %s--------------------------"%step)
            command="%s %s -jar %s BaseRecalibrator -R %s -I %s  -O %s/recal_data.table %s" % (
                java, java_opts, gatk, ref_genome,split_bam,work_gatk,BaseRecalibrator_opts)
            command="bash -c \"%s\""%command      
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
        else:
            logger.info("Skipping step %d: %s"%(step,msg))
        step+=1

        msg = "GATK ApplyBQSR for %s"%sample
        if start<=step:
            logger.info("--------------------------STEP %s--------------------------"%step)
            command="%s %s -jar %s ApplyBQSR -R %s -I %s -bqsr %s/recal_data.table -O %s/bsqr.bam %s" % (
                java, java_opts, gatk, ref_genome,split_bam,work_gatk,work_gatk,ApplyBQSR_opts)
            command="bash -c \"%s\""%command      
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
        else:
            logger.info("Skipping step %d: %s"%(step,msg))
        step+=1
        split_bam="%s/bsqr.bam"%work_gatk
    else:
        msg = "GATK BaseRecalibrator for %s"%sample
        logger.info("Skipping step %d: %s"%(step,msg))
        step+=1
        msg = "GATK ApplyBQSR for %s"%sample
        logger.info("Skipping step %d: %s"%(step,msg))
        step+=1

    msg = "GATK HaplotypeCaller for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s -jar %s HaplotypeCaller -R %s -I %s -O %s/variants.vcf %s" % (
            java, java_opts, gatk, ref_genome,split_bam,work_gatk,HaplotypeCaller_opts)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "GATK VariantFiltration for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s -jar %s VariantFiltration -R %s -V %s/variants.vcf -O %s/variants_filtered.vcf %s" % (
            java, java_opts, gatk, ref_genome,work_gatk,work_gatk,VariantFiltration_opts)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    out_gatk=os.path.join(outdir,"gatk",sample)
    create_dirs([out_gatk])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/variants_filtered.vcf"%work_gatk):
            command = "cp %s/variants_filtered.vcf %s/variants_filtered.vcf"%(
                       work_gatk, out_gatk)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=gatk_log_fd, cmd_log=gatk_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    variants = ""
    if os.path.exists("%s/variants_filtered.vcf"%out_gatk):
        logger.info("GATK was successfull!")
        logger.info("Output variants: %s/variants_filtered.vcf"%out_gatk)
        variants = "%s/variants_filtered.vcf"%out_gatk   
    else:            
        logger.info("GATK failed!")
    return variants

def run_variant(variant_caller="GATK", alignment="",
                  ref_genome="", knownsites="",
                  picard=PICARD, gatk=GATK,                  
                  java=JAVA, java_opts="",
                  CleanSam=False, no_BaseRecalibrator=False,                  
                  AddOrReplaceReadGroups_opts="",  MarkDuplicates_opts="",
                  SplitNCigarReads_opts="",
                  BaseRecalibrator_opts="",
                  ApplyBQSR_opts="",  HaplotypeCaller_opts="",
                  VariantFiltration_opts="",  
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT, ignore_exceptions=False):
    variants=""
    if variant_caller.upper()=="GATK":
        try:
            variants=run_gatk(alignment=alignment,
                  ref_genome=ref_genome, knownsites=knownsites,
                  picard=picard, gatk=gatk,                  
                  java=java, java_opts=java_opts,
                  CleanSam=CleanSam,
                  no_BaseRecalibrator=no_BaseRecalibrator,                 
                  AddOrReplaceReadGroups_opts=AddOrReplaceReadGroups_opts,  
                  MarkDuplicates_opts=MarkDuplicates_opts,
                  SplitNCigarReads_opts=SplitNCigarReads_opts,  
                  BaseRecalibrator_opts=BaseRecalibrator_opts,
                  ApplyBQSR_opts=ApplyBQSR_opts,  HaplotypeCaller_opts=HaplotypeCaller_opts,
                  VariantFiltration_opts=VariantFiltration_opts,  
                  start=start, sample= sample, nthreads=nthreads, 
                  workdir=workdir, outdir=outdir, timeout=timeout)
        except Exception as excp:
            logger.info("GATK failed!")
            logger.error(excp)
            if not ignore_exceptions:
                raise Exception(excp)
    return variants
    
    

