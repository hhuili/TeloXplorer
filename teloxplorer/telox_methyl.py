import os,sys
import logging
import time
import resource
import argparse
from pathlib import Path
from . import utilities as utils
from . import telox_asm
from . import MethylPlotter


def run_telox_methyl(args):

    #### Intermediate files
    tel_bam = os.path.join(args.dorado_align_dir, f'{args.out_prefix}.tel.bam')
    aligned_bam = os.path.join(args.dorado_align_dir, f'{args.out_prefix}.tel.aligned.bam')
    sorted_bam = os.path.join(args.dorado_align_dir, f'{args.out_prefix}.tel.sorted.bam')
    bedmethyl = os.path.join(args.tel_methyl_dir, f'{args.out_prefix}.pileup.bed')

    ## found reference genome telomere-subtelomere boundary
    logging.info(f"\n> STEP1: founding telomere boundary")

    path = Path(args.ref_genome)
    assembly_name = path.name.split('.')[0]

    telox_asm.run(
        assembly_input=args.ref_genome,
        outdir=args.boundary_dir,
        out_prefix=assembly_name,
        args=args
    )

    ## subset telomere reads from modBAM
    logging.info(f"\n> STEP2: extracting and mapping telomere reads")
    cmd = ['samtools','view','-b','-N',args.tel_reads,args.modBAM,'-o',tel_bam]
    utils.run_cmd(cmd)

    ## run dorado align
    cmd = ['dorado','aligner','-Y','--threads',str(args.threads),args.ref_genome,tel_bam,'>',aligned_bam]
    utils.run_cmd(cmd,use_shell=True)

    ## sort bam
    cmd = ['samtools','sort',aligned_bam,'-o',sorted_bam]
    utils.run_cmd(cmd)

    ## index bam
    cmd = ['samtools','index',sorted_bam]
    utils.run_cmd(cmd)

    ## run modkit
    logging.info(f"\n> STEP3: converting modBAM to bedMethyl files (Modkit)")
    cmd = ['modkit','pileup',args.modkit_options,sorted_bam,bedmethyl,'--threads',str(args.threads),'--ref',args.ref_genome,'--preset','traditional']
    utils.run_cmd(cmd,use_shell=True)

    ## plot telomere methylation
    if args.plot:
        logging.info(f"\n> STEP4: plotting methylation levels flanking telomere-subtelomere boundaries")
        boundary = os.path.join(args.boundary_dir, f'{assembly_name}.chromtel_length.tsv')

        MethylPlotter.run_methyl_plot(bedmethyl=bedmethyl,
                                      boundary=boundary,
                                      outdir=args.tel_methyl_dir,
                                      out_prefix=args.out_prefix,
                                      args=args)

def setup_common_args(parser):

    group = parser.add_argument_group("General Options")
    group.add_argument("-l", "--tel-reads",
                        required=True,
                        help="Telomere reads file from telox-length (read IDs in first column)")
    group.add_argument("-b", "--modBAM",
                        required=True,
                        help="Modified base BAM file (modBAM)")
    group.add_argument("-R", "--tel-repeats",
                        required=True,
                        help="Telomere repeat definition file (one per line: arm<TAB>repeat, regex supported)")
    group.add_argument("-r", "--ref-genome",
                        required=True,
                        help="Reference genome file")
    group.add_argument("--modkit-options",
                        default="--filter-threshold C:0.85 --mod-threshold m:0.8 --mod-threshold h:0.9",
                        help="Modkit pileup parameters [default: %(default)s]")
    group.add_argument("--plot",
                        action='store_true',
                        help="Generate telomere methylation plot [default: False]")
    group.add_argument("-o", "--out-prefix",
                       default="plot",
                       help="Output prefix [default: %(default)s].")
    group.add_argument("--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")
    group.add_argument("-t", "--threads",
                        type=int,
                        default=4,
                        help="Number of threads for dorado and modkit [default: %(default)s]")

    return parser

def parse_arguments():

    parser = argparse.ArgumentParser(
        description='Telox-methyl: chromosome-specific telomere methylation analysis'
    )

    parser = setup_common_args(parser)

    # Add module-specific arguments
    parser = telox_asm.setup_asm_args(parser)
    parser = MethylPlotter.setup_plotting_args(parser)

    # Apply preset and default options
    args = parser.parse_args()
    args = telox_asm.apply_presets_and_defaults(args)

    return args


def telox_methyl_cli():
    """
    The main entry point for the 'telox-methyl' command-line tool.
    """

    pipeline_start_time = time.time()

    module = "telox-methyl"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    log_file_path = os.path.join(args.outdir, f"{args.out_prefix}.log")
    utils.setup_logging(log_file_path)

    #### Log startup information
    utils.log_parameter_summary(args, module)

    #### Define output paths and file names
    paths = {
        'boundary_dir': "01_ref_boundary",
        'dorado_align_dir': "02_dorado_align",
        'tel_methyl_dir': "03_telomere_methyl",
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    run_telox_methyl(args)

    logging.info(f"\n**** {module} completed! ****")

    total_elapsed = time.time() - pipeline_start_time
    usage = resource.getrusage(resource.RUSAGE_SELF)
    cpu_time_used = usage.ru_utime + usage.ru_stime

    logging.info(f"CPU time(s) : {cpu_time_used:.1f}")
    logging.info(f"Elapsed time(s) : {total_elapsed:.1f}")

if __name__ == "__main__":

    telox_methyl_cli()

