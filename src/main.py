from collections import defaultdict
import sys
import os

from defaults import *
from run_sr_align import run_sr_align
from run_sr_align import run_sr_align
from run_reconstruct import run_reconstruct
from run_quantify import run_quantify
from run_diff import run_diff
from run_dnv_assemebly import run_dnv_assemebly
from run_lr_correct import run_lr_correct
from run_lr_align import run_lr_align
from run_lr_reconstruct import run_lr_reconstruct
from run_lr_fusion import run_lr_fusion
from run_variant import run_variant
from run_editing import run_editing
from run_fusion import run_fusion
from _version import __version__

from utils import *
import logging


def run_pipeline(args,parser):
    create_dirs([args.workdir, args.outdir,os.path.join(args.workdir,"logs")])
    log_file=os.path.join(args.workdir,"logs","run-%s-sample-%s.log"%(time.strftime("%Y%m%d-%H%M%S"),args.sample))
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT, filename=log_file, filemode="w")
    logFormatter = logging.Formatter(FORMAT)
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    logger.info("Running RNASeqPipeline %s" % __version__)
    logger.info("Command-line %s" % (" ".join(sys.argv)))
    logger.info("Arguments are " + str(args))
    logger.info("Run log will be saved in " + log_file)
    mode = args.mode
    logger.info("Run in mode: " + mode)

#     Simple check for arguments
    if mode=="align":
        if (vars(args)["1"]=="" or vars(args)["2"]=="") and args.U=="" and args.sra=="":
            parser.print_help()
            logger.error("Input sequence file(s) are missing.")
            return os.EX_USAGE
    elif mode=="quantify":
        if (vars(args)["1"]=="" or vars(args)["2"]=="") and args.U=="":
            parser.print_help()
            logger.error("Input sequence file(s) are missing.")
            return os.EX_USAGE
    elif mode=="diff":
        if (not args.quant_files or not args.ref_gtf) and \
           (not args.alignments or (not args.transcripts_gtfs and not args.ref_gtf)):
            parser.print_help()
            logger.error("\n\tYou should either provode {the quantification files and a refrence GTF}, \n\
\tOR {the alignment files and a (reference or assembled) GTF files}.")
            return os.EX_USAGE
    elif mode=="denovo":
        if (vars(args)["1"]=="" or vars(args)["2"]=="") and args.U=="" and args.I=="":
            parser.print_help()
            logger.error("Input sequence file(s) are missing.")
            return os.EX_USAGE
    elif mode=="variant":
        if args.no_BaseRecalibrator==False and args.knownsites=="":
            parser.print_help()
            logger.error("\n\tTo run BaseRecalibrator step, knownsites should provide. \n\
\tIf you don't have knownsites, please use --no_BaseRecalibrator option.")
            return os.EX_USAGE




    if mode=="align":
        if not args.sr_aligner.upper()=="HISAT2":
            logger.error("%s is not supported. \
            \nThe supported short read aligner(s) are: %s."%(args.sr_aligner,SR_ALIGNERS))
            return os.EX_USAGE
        logger.info("Assigned sample ID: %s"%args.sample)
        logger.info("Running align step using %s"%args.sr_aligner)
        run_sr_align(sr_aligner=args.sr_aligner, align_idx=args.align_idx,
                      seq_1=vars(args)["1"], seq_2=vars(args)["2"], seq_u=args.U,
                      seq_sra=args.sra, ref_gtf=args.ref_gtf, 
                      hisat2_opts=args.hisat2_opts, hisat2=args.hisat2, hisat2_sps=args.hisat2_sps, samtools=args.samtools,
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="reconstruct":
        if not args.reconstructor.upper()=="STRINGTIE":
            logger.error("%s is not supported. \
            \n The supported transcriptome reconstructor(s) are: %s."%(args.reconstructor,RECONSTRUCTORS))
            return os.EX_USAGE
        logger.info("Assigned sample ID: %s"%args.sample)
        logger.info("Running reconstruct step using %s"%args.reconstructor)
        run_reconstruct(reconstructor=args.reconstructor, alignment_bam=args.alignment_bam,
                      ref_gtf=args.ref_gtf, 
                      stringtie_opts=args.stringtie_opts, stringtie=args.stringtie,
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="quantify":
        if not args.quantifier.upper()=="SALMON-SMEM":
            logger.error("%s is not supported. \
            \nThe supported treanscriptome reconstructor(s) are: %s."%(args.quantifier, QUANTIFIERS))
            return os.EX_USAGE
        logger.info("Assigned sample ID: %s"%args.sample)
        logger.info("Running quantification step using %s"%args.quantifier)
        run_quantify(quantifier=args.quantifier, quantifier_idx=args.quantifier_idx,
                      seq_1=vars(args)["1"], seq_2=vars(args)["2"], seq_u=args.U,
                      salmon_k=args.salmon_k, libtype=args.libtype,                      
                      salmon_smem_opts=args.salmon_smem_opts, salmon=args.salmon,
                      start=args.start, sample= args.sample, nthreads=args.threads, unzip=args.unzip,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="diff":
        if not args.difftool.upper()=="DESEQ2":
            logger.error("%s is not supported. \
            \nThe supported differential analysis tool(s) are: %s."%(args.difftool,DIFFS))
            return os.EX_USAGE
        logger.info("Running differential analysis step using %s"%args.difftool)
        run_diff(difftool=args.difftool, quant_files=args.quant_files, alignments=args.alignments,
                      transcripts_gtfs=args.transcripts_gtfs,
                      ref_gtf=args.ref_gtf,
                      featureCounts_opts=args.featureCounts_opts, featureCounts=args.featureCounts,
                      stringtie=args.stringtie, stringtie_merge_opts=args.stringtie_merge_opts,                  
                      mincount=args.mincount, alpha=args.alpha, 
                      R=args.R, start=args.start, samples=args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)

    elif mode=="denovo":
        if not args.assembler.upper()=="OASES":
            logger.error("%s is not supported. \
            \nThe supported de novo assembler(s) are: %s."%(args.assembler,DNV_ASSEMBLERS))
            return os.EX_USAGE
        logger.info("Running de novo assembly step using %s"%args.assembler)
        run_dnv_assemebly(assembler=args.assembler, assmebly_hash=args.assmebly_hash,
                      seq_1=vars(args)["1"], seq_2=vars(args)["2"], seq_u=args.U, seq_i=args.I,
                      file_format=args.file_format, read_type=args.read_type, 
                      oases=args.oases, velvetg=args.velvetg, velveth=args.velveth,
                      oases_opts=args.oases_opts, velvetg_opts=args.velvetg_opts, velveth_opts=args.velveth_opts,
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="long_correct":
        if not args.long_corrector.upper()=="LORDEC":
            logger.error("%s is not supported. \
            \nThe supported long read error correction tool(s) are: %s."%(args.long_corrector,LR_CORRECTORS))
            return os.EX_USAGE
        logger.info("Running long read error correction step using %s"%args.long_corrector)
        run_lr_correct(long_corrector=args.long_corrector, kmer=args.kmer,
                      solid=args.kmer,long=args.long, short=args.short, 
                      lordec=args.lordec, lordec_opts=args.lordec_opts,
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="long_align":
        if not args.long_aligner.upper()=="STARLONG":
            logger.error("%s is not supported. \
            \nThe supported long read aligner(s) are: %s."%(args.long_aligner,LR_ALIGNERS))
            return os.EX_USAGE
        logger.info("Running long read alignment step using %s"%args.long_aligner)
        run_lr_align(long_aligner=args.long_aligner,long=args.long,
                      genome_dir=args.genome_dir, ref_gtf=args.ref_gtf,
                      starlong=args.starlong, starlong_opts=args.starlong_opts, 
                      sam2psl=args.sam2psl, samtools=args.samtools,
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="long_reconstruct":
        if not args.long_reconstructor.upper()=="IDP":
            logger.error("%s is not supported. \
            \nThe supported long read transcriptome reconstructor(s) are: %s."%(args.long_reconstructor,
            LR_RECONSTRUCTOR))
            return os.EX_USAGE
        logger.info("Running long read transcriptome reconstruction step using %s"%args.long_reconstructor)
        run_lr_reconstruct(long_reconstructor=args.long_reconstructor,
                      alignment=args.alignment, 
                      short_junction=args.short_junction, long_alignment=args.long_alignment,
                      mode_number=args.mode_number,
                      ref_genome=args.ref_genome, ref_all_gpd=args.ref_all_gpd, ref_gpd=args.ref_gpd,
                      read_length=args.read_length,
                      samtools=args.samtools, idp=args.idp, idp_cfg=args.idp_cfg, 
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="long_fusion":
        if not args.long_fusion_caller.upper()=="IDP-FUSION":
            logger.error("%s is not supported. \
            \nThe supported long read fusion detection tool(s)  are: %s."%(args.long_fusion_caller,
            LR_FUSION))
            return os.EX_USAGE
        logger.info("Running long read fusion detection step using %s"%args.long_fusion_caller)
        run_lr_fusion(long_fusion_caller=args.long_fusion_caller,
                      alignment=args.alignment, 
                      short_junction=args.short_junction, long_alignment=args.long_alignment,
                      short_fasta=args.short_fasta, long_fasta=args.long_fasta, 
                      mode_number=args.mode_number,
                      ref_genome=args.ref_genome, ref_all_gpd=args.ref_all_gpd, ref_gpd=args.ref_gpd,
                      uniqueness_bedgraph=args.uniqueness_bedgraph,
                      genome_bowtie2_idx=args.genome_bowtie2_idx, transcriptome_bowtie2_idx=args.transcriptome_bowtie2_idx,
                      read_length=args.read_length,
                      samtools=args.samtools, idpfusion=args.idpfusion, idpfusion_cfg=args.idpfusion_cfg, 
                      gmap=args.gmap, gmap_idx=args.gmap_idx, star_dir=args.star_dir, bowtie2_dir=args.bowtie2_dir,
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="variant":
        if not args.variant_caller.upper()=="GATK":
            logger.error("%s is not supported. \
            \nThe supported variant caller(s) are: %s."%(args.variant_caller,
            variant_caller))
            return os.EX_USAGE
        logger.info("Running variant calling step using %s"%args.variant_caller)
        run_variant(variant_caller=args.variant_caller,
                      alignment=args.alignment, ref_genome=args.ref_genome, knownsites=args.knownsites,
                      picard=args.picard, gatk=args.gatk, 
                      java=args.java, java_opts=args.java_opts,
                      CleanSam=args.CleanSam, IndelRealignment=args.IndelRealignment,
                      no_BaseRecalibrator=args.no_BaseRecalibrator,
                      AddOrReplaceReadGroups_opts=args.AddOrReplaceReadGroups_opts,
                      MarkDuplicates_opts=args.MarkDuplicates_opts, 
                      SplitNCigarReads_opts=args.SplitNCigarReads_opts, 
                      RealignerTargetCreator_opts=args.RealignerTargetCreator_opts, 
                      IndelRealigner_opts=args.IndelRealigner_opts, 
                      BaseRecalibrator_opts=args.BaseRecalibrator_opts, 
                      PrintReads_opts=args.PrintReads_opts, 
                      HaplotypeCaller_opts=args.HaplotypeCaller_opts, 
                      VariantFiltration_opts=args.VariantFiltration_opts, 
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="editing":
        if not args.editing_caller.upper()=="GIREMI":
            logger.error("%s is not supported. \
            \nThe supported RNA editing caller(s) are: %s."%(args.editing_caller,
            editing_caller))
            return os.EX_USAGE
        logger.info("Running RNA editing calling step using %s"%args.editing_caller)
        run_editing(editing_caller=args.editing_caller,
                      alignment=args.alignment, variant=args.variant, 
                      strand_pos=args.strand_pos, genes_pos=args.genes_pos,
                      ref_genome=args.ref_genome, knownsites=args.knownsites,
                      giremi_dir=args.giremi_dir, htslib_dir=args.htslib_dir, 
                      samtools=args.samtools, gatk=args.gatk, 
                      java=args.java, giremi_opts=args.giremi_opts,java_opts=args.java_opts,
                      VariantAnnotator_opts=args.VariantAnnotator_opts, 
                      start=args.start, sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="fusion":
        if not args.fusion_caller.upper()=="FUSIONCATCHER":
            logger.error("%s is not supported. \
            \nThe supported fusion predictor(s) are: %s."%(args.fusion_caller,
            fusion_caller))
            return os.EX_USAGE
        logger.info("Running Fusion prediction step using %s"%args.fusion_caller)
        run_fusion(fusion_caller=args.fusion_caller,
                      data_dir=args.data_dir, input=args.input, start=args.start, 
                      fusioncatcher=args.fusioncatcher, fusioncatcher_opts=args.fusioncatcher_opts, 
                      sample= args.sample, nthreads=args.threads,
                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout)
    elif mode=="all":
        if not args.sr_aligner.upper()=="HISAT2":
            logger.error("%s is not supported. \
            \nThe supported short read aligner(s) are: %s."%(args.sr_aligner,SR_ALIGNERS))
            return os.EX_USAGE
        if not args.reconstructor.upper()=="STRINGTIE":
            logger.error("%s is not supported. \
            \n The supported transcriptome reconstructor(s) are: %s."%(args.reconstructor,RECONSTRUCTORS))
            return os.EX_USAGE
        if not args.quantifier.upper()=="SALMON-SMEM":
            logger.error("%s is not supported. \
            \nThe supported treanscriptome reconstructor(s) are: %s."%(args.quantifier, QUANTIFIERS))
            return os.EX_USAGE
        if not args.difftool.upper()=="DESEQ2":
            logger.error("%s is not supported. \
            \nThe supported differential analysis tool(s) are: %s."%(args.difftool,DIFFS))
            return os.EX_USAGE
        if not args.assembler.upper()=="OASES":
            logger.error("%s is not supported. \
            \nThe supported de novo assembler(s) are: %s."%(args.assembler,DNV_ASSEMBLERS))
            return os.EX_USAGE
        if not args.long_corrector.upper()=="LORDEC":
            logger.error("%s is not supported. \
            \nThe supported long read error correction tool(s) are: %s."%(args.long_corrector,LR_CORRECTORS))
            return os.EX_USAGE
        if not args.long_aligner.upper()=="STARLONG":
            logger.error("%s is not supported. \
            \nThe supported long read aligner(s) are: %s."%(args.long_aligner,LR_ALIGNERS))
            return os.EX_USAGE
        if not args.long_reconstructor.upper()=="IDP":
            logger.error("%s is not supported. \
            \nThe supported long read transcriptome reconstructor(s) are: %s."%(args.long_reconstructor,
            LR_RECONSTRUCTOR))
            return os.EX_USAGE
        if not args.variant_caller.upper()=="GATK":
            logger.error("%s is not supported. \
            \nThe supported variant caller(s) are: %s."%(args.variant_caller,
            variant_caller))
            return os.EX_USAGE
        if not args.editing_caller.upper()=="GIREMI":
            logger.error("%s is not supported. \
            \nThe supported RNA editing caller(s) are: %s."%(args.editing_caller,
            editing_caller))
            return os.EX_USAGE
        if not args.fusion_caller.upper()=="FUSIONCATCHER":
            logger.error("%s is not supported. \
            \nThe supported fusion predictor(s) are: %s."%(args.fusion_caller,
            fusion_caller))
            return os.EX_USAGE
        if not args.long_fusion_caller.upper()=="IDP-FUSION":
            logger.error("%s is not supported. \
            \nThe supported long read fusion detection tool(s)  are: %s."%(args.long_fusion_caller,
            LR_FUSION))
            return os.EX_USAGE


        do_short = True
        if (vars(args)["1"]=="" or vars(args)["2"]=="") and args.U=="":
            parser.print_help()
            logger.info("Input short-read sequence file(s) are missing. Will skipp short-read steps")
            do_short = False
            
        if (vars(args)["1"]=="" and vars(args)["2"]=="") and (args.U==""):
            parser.print_help()
            logger.error("In pipeline mode, only one input short-read type is possible: paired-end (--1 and --2) or unpaired (--U)")
            return os.EX_USAGE
        
        do_long = args.long != ""
        if not do_long:
            logger.info("Input long-read sequence file(s) are missing. Will skipp long-read steps")


        samples=[[replicate for replicate in sample.split(",")] for sample  in args.sample]
        all_samples=[replicate for sample  in samples for replicate in sample]
        n_samples=sum(map(lambda x:len(x),samples))
        input_sr={}
        if (vars(args)["1"] and vars(args)["2"]):
            logger.info("Inputs are paired-end reads.")
            input_sr["1"] = [j for i in vars(args)["1"] for j in i.split(",")]
            input_sr["2"] = [j for i in vars(args)["2"] for j in i.split(",")]
            if len(input_sr["1"])!=n_samples or len(input_sr["2"])!=n_samples:
                parser.print_help()
                logger.error("Number of short paired-end input sequences does not match number of samples.")
                return os.EX_USAGE
            input_mode="paired"
            input_sr["1"]={all_samples[i]:j for i,j in enumerate(input_sr["1"])}
            input_sr["2"]={all_samples[i]:j for i,j in enumerate(input_sr["2"])}
            input_sr["U"]={all_samples[i]:"" for i,j in enumerate(input_sr["1"])}
        else:
            logger.info("Inputs are unpaired reads.")
            input_sr["U"] = [j for i in args.U for j in i.split(",")]
            if len(input_sr["U"])!=n_samples:
                parser.print_help()
                logger.error("Number of short unpiared input sequences does not match number of samples.")
                return os.EX_USAGE
            input_mode="un-paired"
            input_sr["U"]={all_samples[i]:j for i,j in enumerate(input_sr["U"])}
            input_sr["1"]={all_samples[i]:"" for i,j in enumerate(input_sr["U"])}
            input_sr["2"]={all_samples[i]:"" for i,j in enumerate(input_sr["U"])}
        
        input_lr={}
        if do_long:
            input_lr = [j for i in args.long for j in i.split(",")]
            if len(input_lr)!=n_samples:
                parser.print_help()
                logger.error("Number of long input sequences does not match number of samples.")
                return os.EX_USAGE
            input_lr={all_samples[i]:j for i,j in enumerate(input_lr)}

        alignments_bam={}
        junctions_tab={}
        junctions_bed={}
        transcripts={}
        abundances={}
        quant={}
        diff_af=""
        diff_al=""
        variants={}
        transcripts_dnv={}
        edits={}
        fusions={}
        corrected={}
        alignments_lr={}
        transcripts_lr={}
        abundances_lr={}
        fusions_lr={}
        if do_short:
            for si,sample in enumerate(samples):
                alignments_bam[si]={}
                junctions_tab[si]={}
                junctions_bed[si]={}
                transcripts[si]={}
                abundances[si]={}
                quant[si]={}
                transcripts_dnv[si]={}
                for ri,replicate in enumerate(sample):
                    logger.info("Assigned sample ID for replicate-%d in sample-%d: %s"%(ri+1,si+1,replicate))

                    if "align" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running align step using %s for %s"%(args.sr_aligner,replicate))
                        logger.info("******************************************************************************")
                        alignments_bam[si][replicate],junctions_tab[si][replicate],junctions_bed[si][replicate]=run_sr_align(sr_aligner=args.sr_aligner, 
                                      align_idx=args.align_idx,
                                      seq_1=input_sr["1"][replicate], seq_2=input_sr["2"][replicate], 
                                      seq_u=input_sr["U"][replicate],
                                      seq_sra="", ref_gtf=args.ref_gtf, 
                                      hisat2_opts=args.hisat2_opts, hisat2=args.hisat2, 
                                      hisat2_sps=args.hisat2_sps, samtools=args.samtools,
                                      start=0, sample=replicate, nthreads=args.threads,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout,ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding align step using %s for %s"%(args.sr_aligner,replicate))
                        logger.info("******************************************************************************")
                        alignments_bam[si][replicate],junctions_tab[si][replicate],junctions_bed[si][replicate]=["","",""]
                        
                    if "reconstruct" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running reconstruct step using %s for %s"%(args.reconstructor,replicate))
                        logger.info("******************************************************************************")
                        transcripts[si][replicate],abundances[si][replicate]=run_reconstruct(reconstructor=args.reconstructor,
                                      alignment_bam=alignments_bam[si][replicate],
                                      ref_gtf=args.ref_gtf, 
                                      stringtie_opts=args.stringtie_opts, stringtie=args.stringtie,
                                      start=0, sample=replicate, nthreads=args.threads,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding reconstruct step using %s for %s"%(args.reconstructor,replicate))
                        logger.info("******************************************************************************")
                        transcripts[si][replicate],abundances[si][replicate]=["",""]

                    if "quantify" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running quantification step using %s for %s"%(args.quantifier,replicate))
                        logger.info("******************************************************************************")
                        quant[si][replicate]=run_quantify(quantifier=args.quantifier, quantifier_idx=args.quantifier_idx,
                                      seq_1=input_sr["1"][replicate], seq_2=input_sr["2"][replicate], 
                                      seq_u=input_sr["U"][replicate],
                                      salmon_k=args.salmon_k, libtype=args.libtype,                      
                                      salmon_smem_opts=args.salmon_smem_opts, salmon=args.salmon,
                                      start=0, sample=replicate, nthreads=args.threads, unzip=args.unzip,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding quantification step using %s for %s"%(args.quantifier,replicate))
                        logger.info("******************************************************************************")
                        quant[si][replicate]=""

                    if "denovo" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running de novo assembly step using %s for %s"%(args.assembler,replicate))
                        logger.info("******************************************************************************")
                        transcripts_dnv[si][replicate]=run_dnv_assemebly(assembler=args.assembler, 
                                       assmebly_hash=args.assmebly_hash,
                                      seq_1=input_sr["1"][replicate], seq_2=input_sr["2"][replicate], 
                                      seq_u=input_sr["U"][replicate], seq_i="",
                                      file_format=args.file_format, read_type=args.read_type, 
                                      oases=args.oases, velvetg=args.velvetg, velveth=args.velveth,
                                      oases_opts=args.oases_opts, velvetg_opts=args.velvetg_opts, 
                                      velveth_opts=args.velveth_opts,
                                      start=0, sample= replicate, nthreads=args.threads,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding de novo assembly step using %s for %s"%(args.assembler,replicate))
                        logger.info("******************************************************************************")
                        transcripts_dnv[si][replicate]=""

            if "diff" not in args.exclude:
                logger.info("******************************************************************************")
                logger.info("Running differential analysis step (based on alignment-free quantification results) using %s for %s"%(args.difftool,samples))
                logger.info("******************************************************************************")
                diff_af=run_diff(difftool=args.difftool, quant_files=[",".join([quant[si][replicate] for replicate in sample]) for si,sample in enumerate(samples)],
                             alignments="",
                              transcripts_gtfs="",
                              ref_gtf=args.ref_gtf,
                              featureCounts_opts=args.featureCounts_opts, featureCounts=args.featureCounts,
                              stringtie=args.stringtie, stringtie_merge_opts=args.stringtie_merge_opts,                  
                              mincount=args.mincount, alpha=args.alpha, 
                              R=args.R, start=0, samples=args.sample, nthreads=args.threads,
                              workdir=os.path.join(args.workdir, "diff-quant"), 
                              outdir=os.path.join(args.outdir, "diff-quant"), timeout=args.timeout, ignore_exceptions=True)

                logger.info("******************************************************************************")
                logger.info("Running differential analysis step (based on alignment results) using %s for %s"%(args.difftool,samples))
                logger.info("******************************************************************************")
        #         if use_tgtf
                diff_al=run_diff(difftool=args.difftool, quant_files="",
                             alignments=[",".join([alignments_bam[si][replicate] for replicate in sample]) for si,sample in enumerate(samples)],
                              transcripts_gtfs=[",".join([transcripts[si][replicate] for replicate in sample]) for si,sample in enumerate(samples)],
                              ref_gtf=args.ref_gtf,
                              featureCounts_opts=args.featureCounts_opts, featureCounts=args.featureCounts,
                              stringtie=args.stringtie, stringtie_merge_opts=args.stringtie_merge_opts,                  
                              mincount=args.mincount, alpha=args.alpha, 
                              R=args.R, start=0, samples=args.sample, nthreads=args.threads,
                              workdir=os.path.join(args.workdir, "diff-alignment"), 
                              outdir=os.path.join(args.outdir, "diff-alignment"), timeout=args.timeout, ignore_exceptions=True)
            else:
                logger.info("******************************************************************************")
                logger.info("Excluding differential analysis step (based on alignment-free quantification results) using %s for %s"%(args.difftool,samples))
                logger.info("******************************************************************************")
                diff_af=""
                logger.info("******************************************************************************")
                logger.info("Excluding differential analysis step (based on alignment results) using %s for %s"%(args.difftool,samples))
                logger.info("******************************************************************************")
                diff_al=""


            for si,sample in enumerate(samples):
                variants[si]={}
                edits[si]={}
                fusions[si]={}
                for ri,replicate in enumerate(sample):
                    if "variant" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running variant calling step using %s for %s"%(args.variant_caller,replicate))
                        logger.info("******************************************************************************")
                        variants[si][replicate]=run_variant(variant_caller=args.variant_caller,
                                      alignment=alignments_bam[si][replicate], ref_genome=args.ref_genome, 
                                      knownsites=args.knownsites,
                                      picard=args.picard, gatk=args.gatk, 
                                      java=args.java, java_opts=args.java_opts,
                                      CleanSam=args.CleanSam, IndelRealignment=args.IndelRealignment,
                                      no_BaseRecalibrator=args.no_BaseRecalibrator,
                                      AddOrReplaceReadGroups_opts=args.AddOrReplaceReadGroups_opts,
                                      MarkDuplicates_opts=args.MarkDuplicates_opts, 
                                      SplitNCigarReads_opts=args.SplitNCigarReads_opts, 
                                      RealignerTargetCreator_opts=args.RealignerTargetCreator_opts, 
                                      IndelRealigner_opts=args.IndelRealigner_opts, 
                                      BaseRecalibrator_opts=args.BaseRecalibrator_opts, 
                                      PrintReads_opts=args.PrintReads_opts, 
                                      HaplotypeCaller_opts=args.HaplotypeCaller_opts, 
                                      VariantFiltration_opts=args.VariantFiltration_opts, 
                                      start=0, sample=replicate, nthreads=args.threads,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding variant calling step using %s for %s"%(args.variant_caller,replicate))
                        logger.info("******************************************************************************")
                        variants[si][replicate]=""
                
                    if "editing" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running RNA editing calling step using %s for %s"%(args.editing_caller,replicate))
                        logger.info("******************************************************************************")
                        edits[si][replicate]=run_editing(editing_caller=args.editing_caller,
                                      alignment=alignments_bam[si][replicate], variant=variants[si][replicate], 
                                      strand_pos=args.strand_pos, genes_pos=args.genes_pos,
                                      ref_genome=args.ref_genome, knownsites=args.knownsites,
                                      giremi_dir=args.giremi_dir, htslib_dir=args.htslib_dir, 
                                      samtools=args.samtools, gatk=args.gatk, 
                                      java=args.java, giremi_opts=args.giremi_opts,java_opts=args.java_opts,
                                      VariantAnnotator_opts=args.VariantAnnotator_opts, 
                                      start=0, sample= replicate, nthreads=args.threads,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding RNA editing calling step using %s for %s"%(args.editing_caller,replicate))
                        logger.info("******************************************************************************")
                        edits[si][replicate]=""

                    if "fusion" not in args.exclude:
                        logger.info("******************************************************************************")
                        logger.info("Running Fusion prediction step using %s for %s"%(args.fusion_caller,replicate))
                        logger.info("******************************************************************************")
                        fusions[si][replicate]=run_fusion(fusion_caller=args.fusion_caller,
                                      data_dir=args.data_dir, input="%s,%s"%(input_sr["1"][replicate],
                                      input_sr["2"][replicate]) if input_mode=="paired" else input_sr["U"][replicate], 
                                      start=0, 
                                      fusioncatcher=args.fusioncatcher, fusioncatcher_opts=args.fusioncatcher_opts, 
                                      sample= replicate, nthreads=args.threads,
                                      workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                    else:
                        logger.info("******************************************************************************")
                        logger.info("Excluding RNA editing calling step using %s for %s"%(args.editing_caller,replicate))
                        logger.info("******************************************************************************")
                        fusions[si][replicate]=""

        if do_long:
            if do_short:
                for si,sample in enumerate(samples):
                    corrected[si]={}
                    alignments_lr[si]={}
                    transcripts_lr[si]={}
                    abundances_lr[si]={}
                    fusions_lr[si]={}
                    for ri,replicate in enumerate(sample):
                    
                        if "long_correct" not in args.exclude:
                            logger.info("******************************************************************************")
                            logger.info("Running long read error correction step using %s for %s"%(args.long_corrector,replicate))
                            logger.info("******************************************************************************")
                            corrected[si][replicate]=run_lr_correct(long_corrector=args.long_corrector, kmer=args.kmer,
                                          solid=args.kmer,long=input_lr[replicate], short="%s,%s"%(input_sr["1"][replicate],
                                          input_sr["2"][replicate]) if input_mode=="paired" else input_sr["U"][replicate], 
                                          lordec=args.lordec, lordec_opts=args.lordec_opts,
                                          start=0, sample= replicate, nthreads=args.threads,
                                          workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                        else:
                            logger.info("******************************************************************************")
                            logger.info("Excluding long read error correction step using %s for %s"%(args.long_corrector,replicate))
                            logger.info("******************************************************************************")
                            corrected[si][replicate]=""
                    
                    
                    
                        if "long_align" not in args.exclude:
                            logger.info("******************************************************************************")
                            logger.info("Running long read alignment step on corrected long-reads using %s for %s"%(args.long_aligner,replicate))
                            logger.info("******************************************************************************")
                            alignments_lr[si][replicate]=run_lr_align(long_aligner=args.long_aligner,long=corrected[si][replicate],
                                          genome_dir=args.star_genome_dir, ref_gtf=args.ref_gtf,
                                          starlong=args.starlong, starlong_opts=args.starlong_opts, 
                                          sam2psl=args.sam2psl, samtools=args.samtools,
                                          start=0, sample= replicate, nthreads=args.threads,
                                          workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                        else:
                            logger.info("******************************************************************************")
                            logger.info("Excluding long read alignment step on corrected long-reads using %s for %s"%(args.long_aligner,replicate))
                            logger.info("******************************************************************************")
                            alignments_lr[si][replicate]=""

                        if "long_reconstruct" not in args.exclude:
                            logger.info("******************************************************************************")
                            logger.info("Running long read transcriptome reconstruction step using %s for %s"%(args.long_reconstructor,replicate))
                            logger.info("******************************************************************************")
                            transcripts_lr[si][replicate],abundances_lr[si][replicate]=run_lr_reconstruct(long_reconstructor=args.long_reconstructor,
                                          alignment=alignments_bam[si][replicate], 
                                          short_junction=junctions_bed[si][replicate], 
                                          long_alignment=alignments_lr[si][replicate],
                                          mode_number=args.mode_number,
                                          ref_genome=args.ref_genome, ref_all_gpd=args.ref_all_gpd, ref_gpd=args.ref_gpd,
                                          read_length=args.read_length,
                                          samtools=args.samtools, idp=args.idp, idp_cfg=args.idp_cfg, 
                                          start=0, sample= replicate, nthreads=args.threads,
                                          workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                        else:
                            logger.info("******************************************************************************")
                            logger.info("Excluding long read transcriptome reconstruction step using %s for %s"%(args.long_reconstructor,replicate))
                            logger.info("******************************************************************************")
                            transcripts_lr[si][replicate],abundances_lr[si][replicate]=["",""]
                        
                        if "long_fusion" not in args.exclude:
                            logger.info("******************************************************************************")
                            logger.info("Running long read fusion detection step using %s for %s"%(args.long_reconstructor,replicate))
                            logger.info("******************************************************************************")
                            transcripts_lr[si][replicate],abundances_lr[si][replicate]=run_lr_reconstruct(long_reconstructor=args.long_reconstructor,
                                          alignment=alignments_bam[si][replicate], 
                                          short_junction=junctions_bed[si][replicate], 
                                          long_alignment=alignments_lr[si][replicate],
                                          mode_number=args.mode_number,
                                          ref_genome=args.ref_genome, ref_all_gpd=args.ref_all_gpd, ref_gpd=args.ref_gpd,
                                          read_length=args.read_length,
                                          samtools=args.samtools, idp=args.idp, idp_cfg=args.idp_cfg, 
                                          start=0, sample= replicate, nthreads=args.threads,
                                          workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                            fusions_lr[si][replicate]=run_lr_fusion(long_fusion_caller=args.long_fusion_caller,
                                          alignment=alignments_bam[si][replicate], 
                                          short_junction=junctions_bed[si][replicate], 
                                          short_fasta=input_sr["U"][replicate], long_fasta=corrected[si][replicate], 
                                          mode_number=args.mode_number,
                                          ref_genome=args.ref_genome, ref_all_gpd=args.ref_all_gpd, ref_gpd=args.ref_gpd,
                                          uniqueness_bedgraph=args.uniqueness_bedgraph,
                                          genome_bowtie2_idx=args.genome_bowtie2_idx, transcriptome_bowtie2_idx=args.transcriptome_bowtie2_idx,
                                          read_length=args.read_length,
                                          samtools=args.samtools, idpfusion=args.idpfusion, idpfusion_cfg=args.idpfusion_cfg, 
                                          gmap=args.gmap, gmap_idx=args.gmap_idx, star_dir=args.star_dir, bowtie2_dir=args.bowtie2_dir,
                                          start=0, sample= args.sample, nthreads=args.threads,
                                          workdir=args.workdir, outdir=args.outdir, timeout=args.timeout,ignore_exceptions=True)
                        else:
                            logger.info("******************************************************************************")
                            logger.info("Excluding long read transcriptome reconstruction step using %s for %s"%(args.long_reconstructor,replicate))
                            logger.info("******************************************************************************")
                            transcripts_lr[si][replicate],abundances_lr[si][replicate]=["",""]
            else:
                for si,sample in enumerate(samples):
                    corrected[si]={}
                    for ri,replicate in enumerate(sample):
                        if "long_align" not in args.exclude:
                            logger.info("******************************************************************************")
                            logger.info("Running long read alignment step on original long-reads using %s for %s"%(args.long_aligner,replicate))
                            logger.info("******************************************************************************")
                            alignments_lr[si][replicate]=run_lr_align(long_aligner=args.long_aligner,long=input_lr[replicate],
                                          genome_dir=args.star_genome_dir, ref_gtf=args.ref_gtf,
                                          starlong=args.starlong, starlong_opts=args.starlong_opts, 
                                          sam2psl=args.sam2psl, samtools=args.samtools,
                                          start=0, sample= replicate, nthreads=args.threads,
                                          workdir=args.workdir, outdir=args.outdir, timeout=args.timeout, ignore_exceptions=True)
                        else:
                            logger.info("******************************************************************************")
                            logger.info("Excluding long read alignment step on original long-reads using %s for %s"%(args.long_aligner,replicate))
                            logger.info("******************************************************************************")
                            alignments_lr[si][replicate]=""

                    
                                      
        tasks={"Short-read alignment":[alignments_bam,junctions_tab,junctions_bed],
               "Short-read transcriptome reconstruction":[transcripts,abundances],
               "Short-read alignment-free quantification":[quant],
               "Short-read alignment-free differential analysis":[diff_af],
               "Short-read alignment-based differential analysis":[diff_al],
               "Short-read de novo assembly":[transcripts_dnv],
               "Short-read variant calling":[variants],
               "Short-read rna editing detection":[edits],
               "Short-read fusion detection":[fusions],
               "Long-read error correction":[corrected],
               "Long-read alignment":[alignments_lr],
               "long-read transcriptome reconstruction":[transcripts_lr,abundances_lr],
               "long-read fusion detection":[fusions_lr],
        }
        ordered_tasks=["Short-read alignment",
               "Short-read transcriptome reconstruction",
               "Short-read alignment-free quantification",
               "Short-read alignment-free differential analysis",
               "Short-read alignment-based differential analysis",
               "Short-read de novo assembly",
               "Short-read variant calling",
               "Short-read rna editing detection",
               "Short-read fusion detection",
               "Long-read error correction",
               "Long-read alignment",
               "long-read transcriptome reconstruction",
               "long-read fusion detection"]
        success={task:[] for task in ordered_tasks}
        failure={task:[] for task in ordered_tasks}
        for t,vv in tasks.iteritems():
            v=vv[0]
            if t=="Short-read alignment-free differential analysis" or t=="Short-read alignment-based differential analysis":
                if v:
                    success[t].append("ALL")
                else:
                    failure[t].append("ALL")
            else:
                if v:
                    for si,sample in enumerate(samples):
                        for replicate in sample:
                            if v[si][replicate]:
                                success[t].append(replicate)
                            else:
                                failure[t].append(replicate)
                else:
                    failure[t].append("ALL")
        
            
        
        logger.info("***********************************************")
        logger.info("Successfull Runs:")
        logger.info("***********************************************")
        for t in ordered_tasks:
            if not set(success[t])^set(all_samples):
                success[t]=["ALL"]
            if success[t]:
                logger.info("%s: %s"%(t,",".join(success[t])))
        logger.info("")

        logger.info("***********************************************")
        logger.info("Failed Runs:")
        logger.info("***********************************************")
        for t in ordered_tasks:
            if not set(failure[t])^set(all_samples):
                failure[t]=["ALL"]
            if failure[t]:
                logger.info("%s: %s"%(t,",".join(failure[t])))
        logger.info("")
    else:
        logger.error("wrong mode %s"%(mode))
        return os.EX_USAGE

    logger.info("All Done!")

    return os.EX_OK
