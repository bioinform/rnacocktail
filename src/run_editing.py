import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *
import pysam
import sys
import csv
import pybedtools

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

def filter_multi_chr_alignments(in_file,out_file):
    curren_read=""
    chrms=set([])
    reads=[]
    infile=pysam.AlignmentFile(in_file, "rb")
    outfile=pysam.AlignmentFile(out_file, "wb",template=infile)
    for read in infile:
        if read.qname !=curren_read:
            if curren_read!="":
                if len(chrms)==1:
                    for r in reads:
                        outfile.write(r)                        
            curren_read=read.qname
            chrms=set([read.tid])
            reads=[read]
        else:
            chrms.add(read.tid)                
            reads.append(read)
    if len(chrms)==1:
        for r in reads:
            outfile.write(r)                        
    outfile.close()




def fix_SNV_no(feature):
    return pybedtools.Interval(feature.chrom, feature.start, feature.end, name="SNV",
                               score=feature.score, strand=".",otherfields=[".","."])

def merge_info_SNV(feature):
    pos=round(min(abs(int(feature[9])-feature.start),
                  abs(int(feature[10])-feature.start))/float(int(feature[10])-int(feature[9])+1)*100)
    isin=1 if ( feature.start>= int(feature[9]) and feature.start<=int(feature[10])) else -1
    pos=pos*isin
    name="%s,%s"%(feature[3],feature[11])
    otherfields= [str(pos),feature[12]]
    return pybedtools.Interval(chrom=feature.chrom,start=feature.start,end=feature.end,name=name,
                               score=feature.score,strand=feature[13],otherfields=otherfields)

def find_SNV_strands(strand_pos_bed,genes_pos_bed,input_annotated_vcf,output_annotated_bed):

    final_fwd=pybedtools.BedTool(strand_pos_bed).filter(lambda x:x.strand=="+").sort()
    final_rev=pybedtools.BedTool(strand_pos_bed).filter(lambda x:x.strand=="-").sort()

    vcf_intervals=[]
    with open(input_annotated_vcf, 'rb') as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
        for x in spamreader:
            if x[0][0]=="#":
                continue
            if x[6]!="PASS":
                continue
            if len(x[3])!=1 or len(x[4])!=1:
                continue
        
            gt=x[9].split(":")[0]
            gt=gt.split("|") if "|" in gt else gt.split("/")
            if gt[0]==gt[1]:
                continue
        
            vcf_intervals.append(pybedtools.Interval(x[0], int(x[1])-1, int(x[1]), name="SNV",
                                   score=1 if "DB" in x[7] else 0, strand=".",otherfields=[".","."]))
    SNV=pybedtools.BedTool(vcf_intervals).sort().saveas()




    for w in [0,10,50,100,200,400,800,1000]:
        if w==0:
            SNV_no=SNV
        SNV_fwd=SNV_no.window(final_fwd,w=w).each(merge_info_SNV).sort().groupby(g=[1,2,3],c=[4,5,6,7,8],o="first,first,first,max,min")
        SNV_fwd1=SNV_no.window(final_fwd,w=w,v=True)
        SNV_fwd=SNV_fwd.cat(SNV_fwd1,postmerge=False).sort()

        SNV_rev=SNV_no.window(final_rev,w=w).each(merge_info_SNV).sort().groupby(g=[1,2,3],c=[4,5,6,7,8],o="first,first,first,max,min")
        SNV_rev1=SNV_no.window(final_rev,w=w,v=True)
        SNV_rev=SNV_rev.cat(SNV_rev1,postmerge=False).sort()
        SNV_final=SNV_fwd.cat(SNV_rev,postmerge=False).sort().groupby(g=[1,2,3],c=[4,5,6,7,8],o="collapse,first,collapse,collapse,collapse")
    
        SNV_good_=SNV_final.filter(lambda x:len(set(x[5].split(","))-set("."))==1).sort()
        SNV_no=SNV_final.filter(lambda x:len(set(x[5].split(","))-set("."))==0).each(fix_SNV_no).sort()
        SNV_bad_=SNV_final.filter(lambda x:len(set(x[5].split(","))-set("."))>1).sort()
    
        if w==0:
            SNV_good=SNV_good_
            SNV_bad=SNV_bad_
        else:
            SNV_good=SNV_good.cat(SNV_good_,postmerge=False).sort()
            SNV_no=SNV_no.cat(SNV_bad_,postmerge=False).sort()


    SNV_annotated=[]
    cnt=0
    for i in SNV_good:
        name=list(set(i.name.split(","))-set(["SNV"]))[0]
        strand=list(set(i.strand.split(","))-set(["."]))
        strand=strand[0]
        SNV_annotated.append(pybedtools.Interval(chrom=i.chrom,start=i.start,end=i.end,name=name,
                                                 score=i.score,strand=strand))
    for i in SNV_no:
        SNV_annotated.append(pybedtools.Interval(chrom=i.chrom,start=i.start,end=i.end,name="SNV%d"%cnt,
                                                 score=i.score,strand="."))
        cnt+=1
    SNV_output_annotated_bed=pybedtools.BedTool(SNV_annotated).sort()

    Intes=SNV_output_annotated_bed.window(genes_pos_bed,v=True).each(lambda x:
                                    pybedtools.Interval(x[0],int(x[1]),int(x[2]),"Inte",x[4],"#")).sort()
    Genes=SNV_output_annotated_bed.window(genes_pos_bed,u=True)
    SNV_output_annotated_bed=Intes.cat(Genes,postmerge=False).sort().saveas(output_annotated_bed)



def run_giremi(alignment="", variant="", 
                  strand_pos="", genes_pos="",
                  ref_genome="", knownsites="",
                  giremi_dir="", htslib_dir="",
                  samtools=SAMTOOLS, gatk=GATK,                  
                  java=JAVA, giremi_opts="", java_opts="",
                  VariantAnnotator_opts="",  
                  start=0, sample= "", nthreads=1,
                  workdir=None, outdir=None, timeout=TIMEOUT):


    logger.info("Running RNA editing detection (GIREMI) for %s"%sample)
    if not os.path.exists(alignment):
        logger.error("Aborting!")
        raise Exception("No alignment file %s"%alignment)
    if not os.path.exists(variant):
        logger.error("Aborting!")
        raise Exception("No variant VCF file %s"%variant)
    if not os.path.exists(strand_pos):
        logger.error("Aborting!")
        raise Exception("No strand position BED file %s"%strand_pos)
    if not os.path.exists(genes_pos):
        logger.error("Aborting!")
        raise Exception("No genes position BED file %s"%genes_pos)
    if not os.path.exists(ref_genome):
        logger.error("Aborting!")
        raise Exception("No reference genome FASTA file %s"%ref_genome)
    if not os.path.exists(knownsites):
        logger.error("Aborting!")
        raise Exception("No VCF knownsites file %s"%knownsites)
    if giremi_dir:
        if not os.path.exists(giremi_dir):
            logger.error("Aborting!")
            raise Exception("No GIREMI directory %s"%giremi_dir)

    work_giremi=os.path.join(workdir,"giremi",sample)
    create_dirs([work_giremi])
    
    if nthreads>1:
        if "-nt " not in VariantAnnotator_opts:
            VariantAnnotator_opts += " -nt %d"%nthreads 

    if "-Xms" not in java_opts:
        java_opts += " %s"%JAVA_XMS
    if "-Xmx" not in java_opts:
        java_opts += " %s"%JAVA_XMG
    if "-Djava.io.tmpdir" not in java_opts:
        java_opts += " -Djava.io.tmpdir=%s/javatmp/"%(work_giremi)



    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase GIREMI work directory for %s"%sample
        command="rm -rf %s/*" % (
            work_giremi)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg,timeout=timeout)
    step+=1
    
    giremi_log = os.path.join(work_giremi, "giremi.log")
    giremi_log_fd = open(giremi_log, "w")
    
    
    msg = "Sort BAM by name for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s sort -n -@ %d -T %s/alignments.name_sorted -o %s/alignments.name_sorted.bam %s" % (
            samtools, nthreads, work_giremi, work_giremi, alignment)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=giremi_log_fd, cmd_log=giremi_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1
        

    msg = "Filter alignments mapped to multiple chromosoms for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        logger.info(msg)
        filter_multi_chr_alignments("%s/alignments.name_sorted.bam"%work_giremi,"%s/alignments.chr_unique.bam"%work_giremi)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "Sort BAM by pos for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s sort -@ %d -T %s/alignments.pos_sorted -o %s/alignments.pos_sorted.bam %s/alignments.chr_unique.bam" % (
            samtools, nthreads, work_giremi, work_giremi, work_giremi)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=giremi_log_fd, cmd_log=giremi_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg = "GATK VariantAnnotator for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="%s %s -jar %s -T VariantAnnotator -R %s -V %s -L %s -o %s/annotated.vcf --dbsnp %s %s" % (
            java, java_opts, gatk, ref_genome,variant,variant,work_giremi,knownsites,VariantAnnotator_opts)
        command="bash -c \"%s\""%command      
        cmd = TimedExternalCmd(command, logger, raise_exception=True)
        retcode = cmd.run(cmd_log_fd_out=giremi_log_fd, cmd_log=giremi_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    msg="Find variant strands for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        logger.info(msg)
        find_SNV_strands(strand_pos, genes_pos,  "%s/annotated.vcf"%work_giremi, "%s/SNV_annotated.bed"%work_giremi)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    if htslib_dir:
        if "LD_LIBRARY_PATH" in os.environ:
            os.environ["LD_LIBRARY_PATH"] += ":%s/"%htslib_dir
        else:
            os.environ["LD_LIBRARY_PATH"] = htslib_dir

    if giremi_dir:
        os.environ["PATH"] += ":%s/"%giremi_dir
                
    msg = "Run GIREMI for %s"%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        command="cd %s && %s %s -f %s -l %s/SNV_annotated.bed -o %s/giremi_out.txt %s/alignments.pos_sorted.bam" % (
            giremi_dir,GIREMI, giremi_opts, os.path.abspath(ref_genome), os.path.abspath(work_giremi), os.path.abspath(work_giremi),os.path.abspath(work_giremi))
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(cmd_log_fd_out=giremi_log_fd, cmd_log=giremi_log, msg=msg, timeout=timeout)
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

        
    if os.path.exists("%s/giremi_out.txt"%work_giremi) and not os.path.exists("%s/giremi_out.txt.res"%work_giremi):

        msg="Identify N variants for %s"%sample
        if start<=step:
            logger.info("--------------------------STEP %s--------------------------"%step)
            logger.info(msg)
            with open("%s/giremi_out.txt"%work_giremi) as csv_file_i:
                spamreader = csv.reader(csv_file_i, delimiter='\t', quotechar='|')
                with open("%s/N.bed"%work_giremi, 'wb') as csvfile_o:
                    spamwriter = csv.writer(csvfile_o, delimiter='\t',
                                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
                    for row in spamreader:
                        if (row[5]=="N" or row[8]=="N"):
                            spamwriter.writerow([row[0],int(row[1])-1,row[1]])
        else:
            logger.info("Skipping step %d: %s"%(step,msg))
        step+=1

        cnt=len(pybedtools.BedTool("%s/N.bed"%work_giremi))
        if cnt>0:
            msg="Remove N variants for %s"%sample
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
                logger.info(msg)
                pybedtools.BedTool("%s/SNV_annotated.bed"%work_giremi).intersect(
                "%s/N.bed"%work_giremi,r=True, f=1, v=True).saveas("%s/SNV_annotated_filtered.bed"%work_giremi)
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1
            
            msg = "Rerun GIREMI for %s"%sample
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
                if os.path.exists("%s/SNV_annotated_filtered.bed"%work_giremi):
                    command="cd %s && %s %s -f %s -l %s/SNV_annotated_filtered.bed -o %s/giremi_out.txt %s/alignments.pos_sorted.bam" % (
                        giremi_dir,GIREMI, giremi_opts, os.path.abspath(ref_genome), os.path.abspath(work_giremi), os.path.abspath(work_giremi),os.path.abspath(work_giremi))
                    command="bash -c \"%s\""%command        
                    cmd = TimedExternalCmd(command, logger, raise_exception=False)
                    retcode = cmd.run(cmd_log_fd_out=giremi_log_fd, cmd_log=giremi_log, msg=msg, timeout=timeout)
                else:
                    logger.info("No file %s/SNV_annotated_filtered.bed"%work_giremi)
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1
        else:
            step+=2
    else:
        step+=3

    out_giremi=os.path.join(outdir,"giremi",sample)
    create_dirs([out_giremi])
    msg="Copy predictions to output directory for %s."%sample
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if os.path.exists("%s/giremi_out.txt.res"%work_giremi):
            command = "cp %s/giremi_out.txt.res %s/giremi_out.txt.res"%(
                       work_giremi, out_giremi)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=giremi_log_fd, cmd_log=giremi_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1


    edits = ""
    if os.path.exists("%s/giremi_out.txt.res"%out_giremi):
        logger.info("GIREMI was successfull!")
        logger.info("Output edits: %s/giremi_out.txt.res"%out_giremi)
        edits = "%s/giremi_out.txt.res"%out_giremi   
    else:            
        logger.info("GIREMI failed!")
    return edits

def run_editing(editing_caller="GIREMI", alignment="", variant="", 
                  strand_pos="", genes_pos="",
                  ref_genome="", knownsites="",
                  giremi_dir="", htslib_dir="",
                  samtools=SAMTOOLS, gatk=GATK,                  
                  java=JAVA, giremi_opts="", java_opts="",
                  VariantAnnotator_opts="",  
                  start=0, sample= "", nthreads=1, 
                  workdir=None, outdir=None, timeout=TIMEOUT, ignore_exceptions=False):
    edits=""

    if editing_caller.upper()=="GIREMI":
        try:
            edits=run_giremi(alignment=alignment, variant=variant, 
                      strand_pos=strand_pos, genes_pos=genes_pos,
                      ref_genome=ref_genome, knownsites=knownsites,
                      giremi_dir=giremi_dir, htslib_dir=htslib_dir, 
                      samtools=samtools, gatk=gatk,                  
                      java=java, giremi_opts=giremi_opts, java_opts=java_opts,
                      VariantAnnotator_opts=VariantAnnotator_opts,  
                      start=start, sample= sample, nthreads=nthreads, 
                      workdir=workdir, outdir=outdir, timeout=timeout)
        except Exception as excp:
            logger.info("GIREMI failed!")
            logger.error(excp)
            if not ignore_exceptions:
                raise Exception(excp)

    return edits
    
    

