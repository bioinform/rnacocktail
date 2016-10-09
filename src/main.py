from collections import defaultdict
import sys
import logging
import os

from defaults import *
from utils import *
from run_sr_align import run_sr_align
from run_sr_align import run_sr_align
from run_reconstruct import run_reconstruct
from run_quantify import run_quantify
from run_diff import run_diff
from run_dnv_assemebly import run_dnv_assemebly
from run_lr_correct import run_lr_correct
from run_lr_align import run_lr_align
from run_lr_reconstruct import run_lr_reconstruct
from run_variant import run_variant
from run_editing import run_editing
from run_fusion import run_fusion
from _version import __version__

FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)



def run_pipeline(args,parser):
    logger.info("Running RNASeqPipeline %s" % __version__)
    logger.info("Command-line %s" % (" ".join(sys.argv)))
    logger.info("Arguments are " + str(args))
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




#     Create the directories for working
    create_dirs([args.workdir, args.outdir])
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
                      R=args.R, start=args.start, samples=args.samples, nthreads=args.threads,
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
                      samtools=args.samtools, idp=args.idp, idp_cfg=args.idp_cfg, 
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
    else:
        logger.error("wrong mode %s"%(mode))
        return os.EX_USAGE

    logger.info("All Done!")

    return os.EX_OK
