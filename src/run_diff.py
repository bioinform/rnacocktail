import os
from external_cmd import TimedExternalCmd
from defaults import *
from utils import *
import csv
import glob

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logFormatter = logging.Formatter(FORMAT)
logger = logging.getLogger(__name__)
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

def tx2gene_map(ref_gtf_file,tx2gene_file):
    tx2gene={}
    with open(ref_gtf_file, 'r') as input_f:
        for line in input_f:
            if (line[0] == '#'):
                continue
            fields = line.strip().split()
            transcript_info = {k.split()[0]:k.split()[1] for k in ' '.join(fields[8:]).split(";")[:-1]}
            if "transcript_id" in transcript_info and "gene_id" in transcript_info:
                t=transcript_info["transcript_id"].strip().strip("\"")
                g=transcript_info["gene_id"].strip().strip("\"")
                if t not in tx2gene:
                    tx2gene[t]=g
    with open(tx2gene_file , 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(["TXNAME","GENEID"])
        spamwriter.writerows(tx2gene.items())
    return tx2gene

def fix_quant_file(quant_file,fixed_quant_file):
    cnt=0
    with open(quant_file) as q_f:
        with open(fixed_quant_file,"w'") as fixed_q_f:
            for line in q_f:
                if cnt>0:
                    fields=line.strip().split("\t")
                    fields[0]=fields[0].split("|")[0]
                    line="\t".join(fields)+"\n"
                fixed_q_f.write(line)
                cnt+=1

def run_deseq2(quant_files="", alignments="",
              transcripts_gtfs="", ref_gtf="", 
              featureCounts_opts="", featureCounts=FEATURECOUNTS, 
              stringtie=STRINGTIE, stringtie_merge_opts="",                    
              mincount=DESeq2_MINCNT, alpha=DESeq2_ALPHA, R=R_CMD, 
              start=0, samples=[], nthreads=1,
              workdir=None, outdir=None, timeout=TIMEOUT):

    samples=map(lambda x: x.split(","),samples)
    samples_txt="-".join(map(lambda x:",".join(x),samples))

    logger.info("Running differential analysis (DESeq2) for %s"%samples_txt)

    n_samples = len(samples)
    n_replicates=map(len,samples)
    use_quant=True
    use_refgtf=False
    if quant_files and ref_gtf:
        if len(quant_files) != n_samples:
            logger.error("Aborting!")
            raise Exception("Number of input quantification files does not match the number of samples (%s != %s)"%(
            len(quant_files),n_samples))
        quant_files=map(lambda x: x.split(","),quant_files)
        for i,q in enumerate(quant_files):
            if len(q) != n_replicates[i]:
                logger.error("Aborting!")
                raise Exception("Number of input quantification replicate files does not match the number of replicates in %d%s sample  (%s != %s)"%(
                i+1,"st" if i>0 else "th", len(q),n_replicates[i]))
            for r in q:
                if not os.path.exists(r):
                    logger.error("Aborting!")
                    raise Exception("No qantification file %s"%r)
    elif alignments and (transcripts_gtfs or ref_gtf):
        use_quant=False
        if len(alignments) != n_samples:
            logger.error("Aborting!")
            raise Exception("Number of input alignment files does not match the number of samples (%s != %s)"%(
            len(alignments),n_samples))
            
            
        alignments=map(lambda x: x.split(","),alignments)
        for i,a in enumerate(alignments):
            if len(a) != n_replicates[i]:
                logger.error("Aborting!")
                raise Exception("Number of input alignment replicate files does not match the number of replicates in %d%s sample (%s != %s)"%(
                i+1, "st" if i>0 else "th", len(a),n_replicates[i]))
                
            for r in a:
                if not os.path.exists(r):
                    logger.error("Aborting!")
                    raise Exception("No aligment file %s"%r)
        if transcripts_gtfs:
            transcripts_gtfs=map(lambda x: x.split(","),transcripts_gtfs)
            for i,a in enumerate(transcripts_gtfs):
                if len(a) != n_replicates[i]:
                    logger.error("Aborting!")
                    raise Exception("Number of input gtf files does not match the total number of replicates in %d%s sample (%s != %s)"%(
                    i+1, "st" if i>0 else "th", len(a),n_replicates[i]))
        elif ref_gtf:
            use_refgtf=True

        if ref_gtf:   
            if not os.path.exists(ref_gtf):
                logger.error("Aborting!")
                raise Exception("No reference GTF file %s"%ref_gtf)
    else:
        logger.error("Aborting!")
        raise Exception("Either (quantification files + ref_gtf) or (Alignment files + transcripts_gtfs or ref_gtf) is needed.")

    work_deseq2=os.path.join(workdir,"deseq2",samples_txt)
    create_dirs([work_deseq2])

    step=0
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        msg = "Erase DESeq2 work directory for %s"%samples_txt
        command="rm -rf %s/*" % (
            work_deseq2)
        command="bash -c \"%s\""%command        
        cmd = TimedExternalCmd(command, logger, raise_exception=False)
        retcode = cmd.run(msg=msg, timeout=timeout)
    step+=1

    deseq2_log = os.path.join(work_deseq2, "deseq2.log")
    deseq2_log_fd = open(deseq2_log, "w")


    if use_quant:

        msg="prepare tx2gene for %s."%samples_txt
        if start<=step:
            logger.info("--------------------------STEP %s--------------------------"%step)
            tx2gene_file =  ref_gtf.strip() + "tx2gene.csv"
            if os.path.exists(tx2gene_file):
                logger.info("Will use the precomputed %s as tx2gene.csv for %s"%(tx2gene_file,samples_txt))
            else:
                tx2gene_file =  os.path.join(work_deseq2, "tx2gene.csv")
                logger.info("Will computed %s as tx2gene.csv for %s"%(tx2gene_file,samples_txt))
                tx2gene_map(ref_gtf,tx2gene_file)
        else:
            logger.info("Skipping step %d: %s"%(step,msg))
        step+=1
        
        
        msg="compute gene level abundances for %s."%samples_txt
        if start<=step:
            logger.info("--------------------------STEP %s--------------------------"%step)

            fixed_quant_files=[]
            for i,qs in enumerate(quant_files):
                fixed_qs=[]
                for j,q in enumerate(qs):
                    fixed_q =  os.path.join(work_deseq2, "{}.fixed_quant.sf".format(samples[i][j]))
                    fix_quant_file(q,fixed_q)
                    fixed_qs.append(fixed_q)
                fixed_quant_files.append(fixed_qs)

            command = "%s -e \"library('readr'); library('tximport'); \
                       samples=c(%s); (files <- file.path(c(%s))); names(files) <- samples; \
                       tx2gene <- read.csv(file.path('%s'),sep='\\t'); \
                       txi <- tximport(files, type = 'salmon', tx2gene = tx2gene); \
                       save(txi, file='%s/txi.rda'); \
                       write.table(txi$abundance, file = '%s/txi.abundances',\
                       quote = FALSE, sep='\\t'); \
                       write.table(txi$length, file = '%s/txi.length', quote = FALSE, \
                       sep='\\t'); write.table(txi$counts, file = '%s/txi.counts',\
                       quote = FALSE, sep='\\t');\""%(R, ",".join(map(lambda x: "'%s'"%x,reduce(lambda x,y:x+y,samples)))
                                              ,",".join(map(lambda x: "'%s'"%x,reduce(lambda x,y:x+y,fixed_quant_files)))
                                              ,tx2gene_file, work_deseq2, work_deseq2, work_deseq2, work_deseq2)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)
        else:
            logger.info("Skipping step %d: %s"%(step,msg))
        step+=1

        msg = "DESeq2 for %s"%samples_txt
        if start<=step:
            logger.info("--------------------------STEP %s--------------------------"%step)
            command = "%s -e \"library('DESeq2'); load('%s/txi.rda'); \
                       samples <- c(%s); \
                       condition <- factor(c(%s)); \
                       (colData <- data.frame(row.names=colnames(txi$count), condition));\
                        counts <- round(txi$counts); mode(counts) <- 'integer'; \
                        dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~ condition);\
                         stopifnot(txi$countsFromAbundance %%in%% c('no','scaledTPM','lengthScaledTPM')); \
                         if (txi$countsFromAbundance %%in%% c('scaledTPM','lengthScaledTPM')) \
                         {    message('using just counts from tximport');  } else \
                         {    message('using counts and average transcript lengths from tximport'); \
                         lengths <- txi$length;    dimnames(lengths) <- dimnames(dds);\
                         assays(dds)[['avgTxLength']] <- lengths;  }; \
                         dds <- dds[ rowSums(counts(dds)) >= %d, ]; \
                         dds <- DESeq(dds); \
                         for (i in seq_along(condition)){ \
                         for (j in seq_along(condition)){ \
                         if (i < j){\
                         sample1 <- samples[i]; \
                         sample2 <- samples[j]; \
                         res <- results(dds, contrast=c('condition',sample1,sample2), alpha=%f); \
                         (summary(res)); \
                         res_file= sprintf('%s/deseq2_res_%%s_vs_%%s.tab',sample1,sample2);\
                         write.table(res, file = res_file, \
                         quote = FALSE, sep='\\t'); \
                         } \
                         } \
                         } \
                         save(txi,colData,condition,dds,res, \
                         file='%s/deseq2.rda');\""%(
                       R, work_deseq2, ",".join(map(lambda i:"'sample%d'"%(i),range(len(samples)))),
                       ",".join(map(lambda i:"rep('sample%d', %d)"%(i,n_replicates[i]),range(len(samples)))),
                       mincount, alpha, work_deseq2, work_deseq2)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
        else:
            logger.info("Skipping step %d: %s"%(step,msg))
        step+=1
            
    else:
        if use_refgtf:
            msg = "featureCounts for %s"%samples_txt
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
                command="%s %s -o %s/featureCounts.txt -T %d -a %s -g gene_id %s" % (
                featureCounts, featureCounts_opts, work_deseq2, nthreads, ref_gtf, " ".join(reduce(lambda x,y:x+y,alignments)))
                command="bash -c \"%s\""%command
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   

                command="sed -i -e '2s/.*/Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t%s/' %s/featureCounts.txt" % (
                    "\\t".join(reduce(lambda x,y:x+y,samples)),work_deseq2)
                command="bash -c \"%s\""%command
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1


            msg = "DESeq2 for %s"%samples_txt
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
                command = "%s -e \"library('DESeq2'); countData <- read.table('%s/featureCounts.txt', \
                           header=TRUE, row.names=1);  countData <- countData[ ,6:ncol(countData)]; \
                            countData <- as.matrix(countData); \
                           samples <- c(%s); \
                           condition <- factor(c(%s)); \
                           (colData <- data.frame(row.names=colnames(countData), condition));\
                            dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ condition);\
                             dds <- dds[ rowSums(counts(dds)) >= %d, ]; \
                             dds <- DESeq(dds); \
                             for (i in seq_along(condition)){ \
                             for (j in seq_along(condition)){ \
                             if (i < j){\
                             sample1 <- samples[i]; \
                             sample2 <- samples[j]; \
                             res <- results(dds, contrast=c('condition',sample1,sample2), alpha=%f); \
                             (summary(res)); \
                             res_file= sprintf('%s/deseq2_res_%%s_vs_%%s.tab',sample1,sample2);\
                             write.table(res, file = res_file, \
                             quote = FALSE, sep='\\t'); \
                             } \
                             } \
                             } \
                             save(countData,colData,condition,dds,res, \
                             file='%s/deseq2.rda');\""%(
                           R, work_deseq2, ",".join(map(lambda i:"'sample%d'"%(i),range(len(samples)))),
                           ",".join(map(lambda i:"rep('sample%d', %d)"%(i,n_replicates[i]),range(len(samples)))),
                           mincount, alpha, work_deseq2, work_deseq2)
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1

        else:
        
            msg = "Merge transcripts GTFs for %s"%samples_txt
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
            
                if ref_gtf:
                    stringtie_merge_opts += " -G %s"%ref_gtf
                if "-p " not in stringtie_merge_opts:
                    stringtie_merge_opts += " -p %d"%nthreads 
            
                gtfs_list = open("%s/gtfs_list.txt"%work_deseq2, 'w')
                gtfs_list.write("\n".join(reduce(lambda x,y:x+y,transcripts_gtfs)))
                gtfs_list.close()
            
                command="%s --merge %s -o %s/merged.gtf -v %s/gtfs_list.txt" % (
                stringtie, stringtie_merge_opts, work_deseq2, work_deseq2)
                command="bash -c \"%s\""%command
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1

            msg = "featureCounts for %s"%samples_txt
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
                command="%s %s -o %s/featureCounts.txt -T %d -a %s/merged.gtf -g gene_id %s" % (
                featureCounts, featureCounts_opts, work_deseq2, nthreads, work_deseq2, " ".join(reduce(lambda x,y:x+y,alignments)))
                command="bash -c \"%s\""%command
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
            
                command="sed -i -e '2s/.*/Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t%s/' %s/featureCounts.txt" % (
                    "\\t".join(reduce(lambda x,y:x+y,samples)),work_deseq2)
                command="bash -c \"%s\""%command
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1
            
            msg = "DESeq2 for %s"%samples_txt
            if start<=step:
                logger.info("--------------------------STEP %s--------------------------"%step)
                command = "%s -e \"library('DESeq2'); countData <- read.table('%s/featureCounts.txt', \
                           header=TRUE, row.names=1);  countData <- countData[ ,6:ncol(countData)]; \
                            countData <- as.matrix(countData); \
                           samples <- c(%s); \
                           condition <- factor(c(%s)); \
                           (colData <- data.frame(row.names=colnames(countData), condition));\
                            dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~ condition);\
                             dds <- dds[ rowSums(counts(dds)) >= %d, ]; \
                             dds <- DESeq(dds); \
                             for (i in seq_along(condition)){ \
                             for (j in seq_along(condition)){ \
                             if (i < j){\
                             sample1 <- samples[i]; \
                             sample2 <- samples[j]; \
                             res <- results(dds, contrast=c('condition',sample1,sample2), alpha=%f); \
                             (summary(res)); \
                             res_file= sprintf('%s/deseq2_res_%%s_vs_%%s.tab',sample1,sample2);\
                             write.table(res, file = res_file, \
                             quote = FALSE, sep='\\t'); \
                             } \
                             } \
                             } \
                             save(countData,colData,condition,dds, \
                             file='%s/deseq2.rda');\""%(
                           R, work_deseq2, ",".join(map(lambda i:"'sample%d'"%(i),range(len(samples)))),
                           ",".join(map(lambda i:"rep('sample%d', %d)"%(i,n_replicates[i]),range(len(samples)))),
                           mincount, alpha, work_deseq2, work_deseq2)
                cmd = TimedExternalCmd(command, logger, raise_exception=True)
                retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
            else:
                logger.info("Skipping step %d: %s"%(step,msg))
            step+=1

    out_deseq2=os.path.join(outdir,"deseq2",samples_txt)
    create_dirs([out_deseq2])
    msg="Copy predictions to output directory for %s."%samples_txt
    if start<=step:
        logger.info("--------------------------STEP %s--------------------------"%step)
        if len(glob.glob("%s/deseq2_res*.tab"%work_deseq2))>0:
            command = "cp %s/deseq2_res*.tab %s/"%(
                       work_deseq2, out_deseq2)
            cmd = TimedExternalCmd(command, logger, raise_exception=True)
            retcode = cmd.run(cmd_log_fd_out=deseq2_log_fd, cmd_log=deseq2_log, msg=msg, timeout=timeout)   
    else:
        logger.info("Skipping step %d: %s"%(step,msg))
    step+=1

    diff = ""
    if len(glob.glob("%s/deseq2_res*.tab"%out_deseq2))>0:
        logger.info("DESeq2 was successfull!")
        logger.info("Output differential expressions: %s"%(glob.glob("%s/deseq2_res*.tab"%out_deseq2)))
        diff = glob.glob("%s/deseq2_res*.tab"%out_deseq2)  
    else:            
        logger.info("DESeq2 failed!")
    return diff

def run_diff(difftool="DESeq2", quant_files="", alignments="",
              transcripts_gtfs="",
              ref_gtf="",                       
              featureCounts_opts="", featureCounts=FEATURECOUNTS,
              stringtie=STRINGTIE, stringtie_merge_opts="",                    
              mincount=DESeq2_MINCNT, alpha=DESeq2_ALPHA, 
              R=R_CMD, 
              start=0, samples="", nthreads=1,
                  workdir=None, outdir=None, timeout=TIMEOUT, ignore_exceptions=False):
    diff=""
    if difftool.upper()=="DESEQ2":
        try:
            diff=run_deseq2(quant_files=quant_files, alignments=alignments,
                          transcripts_gtfs=transcripts_gtfs, ref_gtf=ref_gtf, 
                          featureCounts_opts=featureCounts_opts, featureCounts=featureCounts,                     
                          stringtie=stringtie, stringtie_merge_opts=stringtie_merge_opts,
                          mincount=mincount, alpha=alpha, 
                          R=R, 
                          start=start, samples=samples, nthreads=nthreads,
                          workdir=workdir, outdir=outdir, timeout=timeout)
        except Exception as excp:
            logger.info("DESeq2 failed!")
            logger.error(excp)
            if not ignore_exceptions:
                raise Exception(excp)
        
    return diff