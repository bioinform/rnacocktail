import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *
import csv

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)



def sort_gpd(in_file,out_file,order_chrs=dict([("%s"%k,k) for k in range(1,23)]+[("MT",23),("X",24),("Y",25)]+[
                                    ("chr%s"%k,k) for k in range(1,23)]+[("chrM",23),("chrX",24),("chrY",25)])):
    with open(in_file) as csv_file:
        spamreader = csv.reader(csv_file, delimiter='\t', quotechar='|')
        rows=[]
        for row in spamreader:
            rows.append(row)
        others_chrs=sorted(set(map(lambda x:x[2],rows))-set(order_chrs.keys()))
        if others_chrs:
            max_id=max(order_chrs.values())
            for i,c in enumerated(others_chrs):
                order_chrs[c]=max_id+i+1
        sorted_rows=sorted(rows,key=lambda x: (order_chrs[x[2]],int(x[4])))
        with open(out_file, 'wb') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter='\t',
                                    quotechar='|', quoting=csv.QUOTE_MINIMAL)
            spamwriter.writerows(sorted_rows)



def run_idp(alignment="", short_junction="", long_alignment="",mode_number=0, 
                  ref_genome="", ref_all_gpd="", ref_gpd="",read_length=100,
                  idp_cfg="", idp=IDP, samtools=SAMTOOLS,
                  start=0, sample= "", nthreads=1,
                  workdir=None, outdir=None, timeout=TIMEOUT):

    logger.info("Running transcriptome lr_reconstruction (IDP) for %s"%sample)
    if not os.path.exists(alignment):
        logger.error("Aborting!")
        raise Exception("No input short read alignment BAM/SAM file %s"%alignment)
    if not os.path.exists(short_junction):
        logger.error("Aborting!")
        raise Exception("No input short read junction BED file %s"%short_junction)
    if not os.path.exists(long_alignment):
        logger.error("Aborting!")
        raise Exception("No input long read alignment PSL file %s"%long_alignment)
        
    if idp_cfg:
        if not os.path.exists(idp_cfg):
            logger.error("Aborting!")
            raise Exception("No input .cfg file %s"%idp_cfg)
        

    
    if mode_number>0:
        start=4
    
    work_idp="%s/idp/%s/"%(workdir,sample)
    create_dirs([work_idp])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase IDP work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_idp)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1



    idp_log = os.path.join(work_idp, "idp.log")
    idp_log_fd = open(idp_log, "w")

    msg = "converting BAM to SAM for %s"%sample
    logger.info("--------------------------STEP %s--------------------------"%step)
    if start<=step:
        if alignment.endswith('.bam'):
            command = "%s view -h -o %s/alignments.sam %s " % (samtools,work_idp,alignment)
            command="bash -c \"%s\""%command       
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=idp_log_fd, cmd_log=idp_log, msg=msg, timeout=timeout)
            alignment =  "%s/alignments.sam"%(work_idp)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    msg = "Preparing run.cfg for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if idp_cfg:
            msg = "copy IDP .cfg file for %s"%sample
            command="cp  %s %s/run.cfg" % (
                idp_cfg, work_idp)
            command="bash -c \"%s\""%command
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=idp_log_fd, cmd_log=idp_log, msg=msg, timeout=timeout)   
        else:
            f=open("%s/run.cfg"%work_idp, 'w')
            f.close()

        cgf_dict={}
        with open("%s/run.cfg"%work_idp , 'r') as cfg_file:
            for line in cfg_file:
                line = line.strip()
                if line=='':
                    continue
                if "=" in line and not line[0]=='#' :
                    k,v=line.split("=")
                    k=k.strip()
                    v=v.strip()
                    cgf_dict[k]=v
                    
        with open("%s/run.cfg"%work_idp , 'w') as cfg_file:
            for k,v in cgf_dict.iteritems():
                cfg_file.write("%s = %s \n"%(k,v))
            if "temp_foldername" not in cgf_dict:
                cfg_file.write("temp_foldername = %s/tmp/ \n"%work_idp)
            if "output_foldername" not in cgf_dict:
                cfg_file.write("output_foldername = %s/out/ \n"%work_idp)
            if "Nthread" not in cgf_dict:
                cfg_file.write("Nthread = %d \n"%nthreads)
            if "LR_psl_pathfilename" not in cgf_dict:
                cfg_file.write("LR_psl_pathfilename = %s \n"%long_alignment)
            if "SR_sam_pathfilename" not in cgf_dict:
                cfg_file.write("SR_sam_pathfilename = %s \n"%alignment)
            if "SR_jun_pathfilename" not in cgf_dict:
                cfg_file.write("SR_jun_pathfilename = %s \n"%short_junction)
            if "genome_pathfilename" not in cgf_dict:       
                cfg_file.write("genome_pathfilename = %s \n"%ref_genome)
            if "allref_annotation_pathfilename" not in cgf_dict:       
                cfg_file.write("allref_annotation_pathfilename = %s \n"%ref_all_gpd)
            if "ref_annotation_pathfilename" not in cgf_dict:       
                cfg_file.write("ref_annotation_pathfilename = %s \n"%ref_gpd)
            if "estimator_choice" not in cgf_dict:       
                cfg_file.write("estimator_choice = MLE \n")
            if "FPR" not in cgf_dict:       
                cfg_file.write("FPR = 0.05 \n")
            if "Njun_limit" not in cgf_dict:       
                cfg_file.write("Njun_limit = 10 \n")
            if "Niso_limit" not in cgf_dict:       
                cfg_file.write("Niso_limit = 100 \n")
            if "aligner_choice" not in cgf_dict:       
                cfg_file.write("aligner_choice = gmap \n")
            if "exon_construction_junction_span" not in cgf_dict:
                cfg_file.write("exon_construction_junction_span = 1 \n")
            if "read_length" not in cgf_dict:
                cfg_file.write("read_length = %d \n"%read_length)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    
    msg = "IDP for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s/run.cfg %d" % (
            idp, work_idp, mode_number)
        command="bash -c \"%s\""%command
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=idp_log_fd, cmd_log=idp_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
    
    msg = "Convert transcript GPD file to GTF for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/out/isoform.gpd"%work_idp):
            sort_gpd("%s/out/isoform.gpd"%work_idp,"%s/out/isoform_sorted.gpd"%work_idp)
            command="gpd2gtf.py \
                  %s/out/isoform_sorted.gpd %s/out/isoform.exp %s/out/isoform.gtf IDP"%(work_idp,work_idp,work_idp)
            command="bash -c \"%s\""%command
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=idp_log_fd, cmd_log=idp_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    out_idp=os.path.join(outdir,"idp",sample)
    create_dirs([out_idp])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/out/isoform.gtf"%work_idp) and \
           os.path.exists("%s/out/isoform.exp"%work_idp):
            command = "cp %s/out/isoform.gtf %s/isoform.gtf"%(
                       work_idp, out_idp)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=idp_log_fd, cmd_log=idp_log, msg=msg, timeout=timeout)   
            
            command = "cp %s/out/isoform.exp %s/isoform.exp"%(
                       work_idp, out_idp)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=idp_log_fd, cmd_log=idp_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1



    transcripts = ""
    abundances = ""
    if os.path.exists("%s/isoform.gtf"%out_idp) and \
       os.path.exists("%s/isoform.exp"%out_idp):
        logger.info("IDP was successfull!")
        logger.info("Output isoforms: %s/isoform.gtf"%out_idp)
        logger.info("Output expressions: %s/isoform.exp"%out_idp)
        transcripts = "%s/isoform.gtf"%out_idp   
        abundances = "%s/isoform.exp"%out_idp   
    else:            
        logger.info("IDP was not successfull!")
    return transcripts,abundances

def run_lr_reconstruct(long_reconstructor="IDP", alignment="",
                  short_junction="", long_alignment="", mode_number=0,
                  ref_genome="", ref_all_gpd="", ref_gpd="", read_length=100,
                  idp_cfg="", idp=IDP, samtools=SAMTOOLS,
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT):
    transcripts = ""
    abundances = ""
    if long_reconstructor.upper()=="IDP":
        transcripts,abundances=run_idp(alignment=alignment, 
                      short_junction=short_junction, long_alignment=long_alignment, 
                      mode_number=mode_number,
                      ref_genome=ref_genome, ref_all_gpd=ref_all_gpd, ref_gpd=ref_gpd,
                      read_length=read_length,
                      idp_cfg=idp_cfg, idp=idp, samtools=samtools,
                      start=start, sample= sample, nthreads=nthreads,
                      workdir=workdir, outdir=outdir, timeout=timeout)
    return transcripts,abundances