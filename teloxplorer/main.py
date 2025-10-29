import os,sys
import logging
import argparse
import time
import resource
from . import telox_align
from . import telox_length
from . import telox_variants
from . import utilities as utils

def setup_common_args(parser):
    
    group = parser.add_argument_group("General Options")
    group.add_argument("-fq", "--long-read-fastq",
                        required=True,
                        help="Input long-read FASTQ file (gzipped or plain)")
    group.add_argument("-r", "--ref-genome",
                        required=True,
                        help="Reference genome file")
    group.add_argument("-R", "--tel-repeats",
                        required=True,
                        help="Telomere repeat definition file (one per line: arm<TAB>repeat, regex supported)")
    group.add_argument("--out-prefix",
                       default="telox",
                       help="Output prefix [default: %(default)s].")
    group.add_argument("--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")
    group.add_argument("-t", "--threads",
                       type=int,
                       default=4,
                       help="Number of threads for minimap2 [default: %(default)s]")
    group.add_argument("--find-TVR",
                       choices=["yes","no"],
                       default="yes",
                       help="Identify telomeric repeat units [default: %(default)s]")
    group.add_argument("--start-step",
                        type=int, 
                        default=1, 
                        choices=[1, 2, 3],
                        help="Step1: extracting telomere-like reads; "
                             "Step2: mapping telomere-like reads; "
                             "Step3: telomere length estimation [default: %(default)s]")
    group.add_argument("--debug",
                       action="store_true",
                       help="Enable debug output")

    return parser

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Teloxplorer: chromosome-specific telomere analysis for long-read sequencing data"
    )

    parser.add_argument("-v", "--version",
                        action="version",
                        version="Teloxplorer v0.0.1",
                        help="Show version and exit")

    parser = setup_common_args(parser)

    # Add module-specific arguments
    parser = telox_align.setup_align_args(parser)
    parser = telox_length.setup_length_args(parser)
    parser = telox_variants.setup_tvr_args(parser)

    # Apply preset and default options
    args = parser.parse_args()
    args = telox_length.apply_presets_and_defaults(args)

    return args

def telox_main_cli():

    pipeline_start_time = time.time()

    module = "teloxplorer"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    log_file_path = os.path.join(args.outdir, f"{args.out_prefix}.log")

    utils.setup_logging(log_file_path, args.debug)

    #### Log startup information
    utils.log_parameter_summary(args, module)

    #### Define output paths and file names
    paths = {
        'tel_reads_dir': "01_telomere_reads",
        'tel_mapping_dir': "02_telomere_mapping",
        'tel_length_dir': "03_telomere_length",
        'tel_variants_dir': "04_telomere_variants"
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    tel_like_fastq = os.path.join(args.tel_reads_dir,f'{args.out_prefix}.tel_like.fastq')

    # Step 1: extracting telomere-like reads
    if args.start_step <= 1:
        logging.info(f"\n> STEP1: extracting telomere-like reads")
        telox_align.extract_telomeric_reads(
            fastq_input=args.long_read_fastq,
            tel_repeats=args.tel_repeats,
            fastq_output=tel_like_fastq
        )
    else:
        logging.info(f"\n> Skipping STEP1")

    # Step 2: alignment-based telomere length estimation
    if args.start_step <= 2:
        logging.info(f"\n> STEP2: mapping telomere-like reads")
        telox_align.run_align(
            fastq_input=tel_like_fastq,
            outdir=args.tel_mapping_dir,
            out_prefix=args.out_prefix,
            args=args
        )
    else:
        logging.info(f"\n> Skipping STEP2")

    # Step 3: alignment-based telomere length estimation
    if args.start_step <= 3:
        logging.info(f"\n> STEP3: alignment-based telomere length estimation")
        bam_file = os.path.join(args.tel_mapping_dir, f"{args.out_prefix}.sort.bam")
        telox_length.run(
            bam_file=bam_file,
            outdir=args.tel_length_dir,
            out_prefix=args.out_prefix,
            args=args
        )
    else:
        logging.info(f"\n> Skipping STEP3")

    # Step 4: Telomere variant repeats analysis
    if args.find_TVR == "yes":
        logging.info(f"\n> STEP4: telomere variant repeats analysis")
        for tel_group in ["chromtel","neotel"]:

            tvr_out_prefix = f"{args.out_prefix}.{tel_group}"
            telomere_fasta = os.path.join(args.tel_length_dir, f"{tvr_out_prefix}_seq.fasta")

            logging.info(f"> Finding TVR for {tel_group}...")
            telox_variants.run(
                fasta_input=telomere_fasta,
                outdir=args.tel_variants_dir,
                out_prefix=tvr_out_prefix,
                args=args)

    logging.info(f"\n**** {module} completed! ****")
    total_elapsed = time.time() - pipeline_start_time

    usage = resource.getrusage(resource.RUSAGE_SELF)
    cpu_time_used = usage.ru_utime + usage.ru_stime
    logging.info(f"CPU time(s) : {cpu_time_used:.1f}")
    logging.info(f"Elapsed time(s) : {total_elapsed:.1f}")

if __name__ == "__main__":

    telox_main_cli()

