import os,sys
import logging
import argparse
from . import utilities as utils
from .telox_align import extract_telomeric_reads,run_telox_align
from .telox_length import run_telox_length
from .telox_variants import run_telox_variants

PRESETS = {
    'Human': {
        'min_tandem_repeats': 5,
        'chrom_subtelomere_boundary': 500000,
        'min_read_length': 1000,
        'min_mapping_quality': 20,
        'min_telomere_freq': 0.6,
        'max_repeat_mismatches': 2,
        'min_initial_telomere_offset': 200,
        'min_telomere_length': 100,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'kmers': "5,6,7",
        'min_variant_repeats': 2,
        'primary_merge': False,
    },
    'Mouse': {
        'min_tandem_repeats': 5,
        'chrom_subtelomere_boundary': 500000,
        'min_read_length': 1000,
        'min_mapping_quality': 20,
        'min_telomere_freq': 0.6,
        'max_repeat_mismatches': 2,
        'min_initial_telomere_offset': 200,
        'min_telomere_length': 100,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'kmers': "5,6,7",
        'min_variant_repeats': 2,
        'primary_merge': False,
    },
    'Yeast': {
        'min_tandem_repeats': 5,
        'chrom_subtelomere_boundary': 20000,
        'min_read_length': 1000,
        'min_mapping_quality': 20,
        'min_telomere_freq': 0.7,
        'max_repeat_mismatches': 0,
        'min_initial_telomere_offset': 100,
        'min_telomere_length': 30,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 50 --ydis 0.2",
        'kmers': "3,4",
        'min_variant_repeats': 4,
        'primary_merge': True,
    }
}

def parse_arguments():
    """Parse command-line arguments for the Teloxplorer pipeline.

    Returns:
        Parsed arguments as a Namespace object.
    """
    parser = argparse.ArgumentParser(
        description="Teloxplorer pipeline for chromosome-specific telomere analysis from long-read sequencing data."
    )
    parser.add_argument(
        "--preset",
        type=str,
        choices=['Human', 'Mouse', 'Yeast'],
        required=False,
        help="Use a preset configuration for a specific species (e.g., Human, Mouse, Yeast). "
            "Manually specified parameters will override preset values."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads for mapping (default: 4).",
    )
    parser.add_argument(
        "-v", "--version",
        action="version",
        version="Teloxplorer v1.0.0",
        help="Show the version number and exit."
    )

    # Group 1: Input Files
    input_group = parser.add_argument_group("Input Files", "Input files for the pipeline.")
    input_group.add_argument(
        "-i", "--long-read-fastq",
        required=True,
        help="Path to the long-read FASTQ file (gzipped or plain).",
    )
    input_group.add_argument(
        "-g", "--ref-genome",
        required=True,
        help="Path to the reference genome file.",
    )
    input_group.add_argument(
        "-r", "--tel-repeats",
        required=True,
        help="Path to telomere repeat definition file (one repeat per line, format: arm<TAB>sequence, supports regex).",
    )

    # Group 2: Telomere Align Parameters
    align_group = parser.add_argument_group("Telomere Align Parameters", "Settings for telomeric read alignment with Minimap2.")
    align_group.add_argument(
        "--min-tandem-repeats",
        type=int,
        default=5,
        help="Minimum consecutive repeat units to identify telomeric reads (default: 5).",
    )
    align_group.add_argument(
        "--minimap2-preset",
        default="-ax map-ont",
        choices=["-ax map-ont", "-ax map-pb", "-ax map-hifi"],
        help="Minimap2 preset option (default: -ax map-ont).",
    )
    align_group.add_argument(
        "-q", "--min-mapping-quality",
        type=int,
        default=20,
        help="Minimum mapping quality (default: 20).",
    )

    # Group 3: Telomere Length Parameters
    telomere_group = parser.add_argument_group("Telomere Length Parameters", "Parameters for calculating telomere length.")
    telomere_group.add_argument(
        "--chrom-subtelomere-boundary",
        type=int,
        default=500000,
        help="Telomere-subtelomere boundary of chromosomes (default: 500000 bp).",
    )
    telomere_group.add_argument(
        "--min-read-length",
        type=int,
        default=1000,
        help="Minimum read length for telomere analysis (default: 1000 bp).",
    )
    telomere_group.add_argument(
        "--min-telomere-freq",
        type=float,
        default=0.6,
        help="Initial scan: frequency threshold for distinguishing telomere from non-telomere segments (default: 0.6).",
    )
    telomere_group.add_argument(
        "--primary-merge",
        action='store_true',
        required=False,
        help="Merge adjacent telomere segments before running BLOOM. (Default: False).",
    )
    telomere_group.add_argument(
        "--max-repeat-mismatches",
        type=int,
        default=2,
        help="Secondary scan: mismatch tolerance for re-labeling non-telomere segments within telomere segments (default: 2).",
    )
    telomere_group.add_argument(
        "--min-initial-telomere-offset",
        type=int,
        default=200,
        help="Minimum non-telomere sequence at the head of a telomere read (default: 200 bp).",
    )
    telomere_group.add_argument(
        "--min-telomere-length",
        type=int,
        default=100,
        help="Minimum telomere length for a telomere read (default: 100 bp).",
    )
    telomere_group.add_argument(
        "--min-subtelomere-length",
        type=int,
        default=100,
        help="Minimum sub-telomere length for a telomere read (default: 100 bp).",
    )

    # Group 4: BLOOM Parameters
    bloom_group = parser.add_argument_group("BLOOM Parameters", "Options for the BLOOM tool.")
    bloom_group.add_argument(
        "--bloom-options",
        type=str,
        default="--bloom_pwidth 500 --bloom_ydis 0.4",
        help="BLOOM preset options (default: --bloom_pwidth 500 --bloom_ydis 0.4)."
            "Used for telomere segments merge.",
    )

    # Group 5: Telomere Variant Repeats Parameters
    variants_group = parser.add_argument_group("Telomere Variant Repeats Parameters", "Parameters for telomere variant repeats idetification.")
    variants_group.add_argument(
        "-k", "--kmers",
        type=str,
        default="5,6,7",
        help="Kmer sizes to find telomeric repeat units (default: 5,6,7).",
    )
    variants_group.add_argument(
        "--min-variant-repeats",
        type=int,
        default=2,
        help="Minimum number of consecutive repeats for telomere variant repeats calling (default: 2).",
    )
    variants_group.add_argument(
        "--find-singleton",
        action='store_true',
        required=False,
        help="Find singleton telomeric repeat units.",
    )

    # Group 6: Output and Debugging
    output_group = parser.add_argument_group("Output and Debugging", "Output directory and debug settings.")
    output_group.add_argument(
        "-o", "--out-prefix",
        default="telox",
        help="Output prefix (default: telox).",
    )
    output_group.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: .).",
    )
    output_group.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug mode for verbose logging.",
    )
    output_group.add_argument(
        '--start-step', 
        type=int, 
        default=1, 
        choices=[1, 2, 3, 4],
        help='Starting step for the pipeline. Defaults to 1 (run all steps).',
    )

    args = parser.parse_args()
    user_provided_args = set()
    for action in parser._actions:
        if action.dest != 'help' and hasattr(args, action.dest):
            if getattr(args, action.dest) != action.default:
                user_provided_args.add(action.dest)

    args.user_provided_args = user_provided_args

    return args

def run_teloxplorer():

    utils.setup_logging()
    args = parse_arguments()

    if args.preset:
        logging.info(f"Applying preset for '{args.preset}':")
        preset_values = PRESETS.get(args.preset)
        if not preset_values:
            available_presets = list(PRESETS.keys())
            raise ValueError(f"Invalid preset '{args.preset}'. Available presets are: {available_presets}")

        for key, value in preset_values.items():
            if key not in args.user_provided_args and hasattr(args, key):
                setattr(args, key, value)
                logging.info(f"  -> Setting '{key}' to preset value: {value}", extra={'notime': True})


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

    #### Intermediate files
    tel_reads_fastq = os.path.join(args.tel_reads_dir, f'{args.out_prefix}.tel_like.fq')
    aln_file = os.path.join(args.tel_mapping_dir, f'{args.out_prefix}.sort.bam')

    #### Run Teloxplorer Pipeline
    logging.info(f"--- Starting pipeline from Step {args.start_step} ---")

    # Step 1: Extracting telomere-like reads
    if args.start_step <= 1:
        logging.info(f"Step 1: Extracting telomere-like reads")
        extract_telomeric_reads(input_fastq=args.long_read_fastq,
                                tel_repeats=args.tel_repeats,
                                min_tandem_repeats=args.min_tandem_repeats,
                                output_fastq=tel_reads_fastq)
    else:
        logging.info(f"Skipping Step 1.")

    # Step 2: Mapping telomere-like reads
    if args.start_step <= 2:
        logging.info(f"Step 2: Mapping telomere-like reads")
        run_telox_align(input_fastq=tel_reads_fastq,
                        ref_genome=args.ref_genome,
                        mm2_preset=args.minimap2_preset,
                        min_mapping_quality=args.min_mapping_quality,
                        out_prefix=args.out_prefix,
                        output_dir=args.tel_mapping_dir,
                        threads=args.threads)
    else:
        logging.info(f"Skipping Step 2.")

    # Step 3: Calculating telomere length
    if args.start_step <= 3:
        logging.info(f"Step 3: Calculating telomere length")
        run_telox_length(aln_file=aln_file,
                        tel_repeats=args.tel_repeats,
                        ref_genome=args.ref_genome,
                        bloom_options=args.bloom_options,
                        out_prefix=args.out_prefix,
                        outdir=args.tel_length_dir,
                        args=args)
    else:
        logging.info(f"Skipping Step 3.")

    # Step 4: Analyzing telomere variant repeats
    if args.start_step <= 4:
        logging.info(f"Step 4: Analyzing telomere variant repeats")
        for tel_group in ['chromtel','neotel']:
            logging.info(f"Finding TVR for {tel_group}...")
            variants_out_prefix = f'{args.out_prefix}.{tel_group}'
            telomere_fasta_file = os.path.join(args.tel_length_dir, f"{variants_out_prefix}.tel_seq.fa")
            run_telox_variants(telomere_fasta=telomere_fasta_file,
                            tel_repeats=args.tel_repeats,
                            out_prefix=variants_out_prefix,
                            outdir=args.tel_variants_dir,
                            args=args)
    else:
        logging.info(f"Skipping Step 4.")

    logging.info(f"--- Teloxplorer pipeline completed successfully! ---\n")

def main_cli():
    """
    The main entry point for the 'teloxplorer' command-line tool.
    """

    run_teloxplorer()

if __name__ == "__main__":

    main_cli()

