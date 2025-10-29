import os
import logging
import re
import pyfastx
import argparse
from . import utilities as utils
from .TelomereMapper import TelomereMapper

def extract_telomeric_reads(fastq_input, tel_repeats, fastq_output, min_tel_repeats=5):

    tel_repeat_dict = utils.parse_tel_repeats(tel_repeats)
    search_str_L = tel_repeat_dict["L"] * min_tel_repeats
    search_str_R = tel_repeat_dict["R"] * min_tel_repeats

    if utils.is_dna_sequence(tel_repeat_dict["L"]):
        search_func = lambda seq: search_str_L in seq or search_str_R in seq
    else:
        pattern_L = re.compile(search_str_L)
        pattern_R = re.compile(search_str_R)
        search_func = lambda seq: pattern_L.search(seq) or pattern_R.search(seq)

    total_reads_found = 0
    batch_size = 1000
    reads_batch = []

    fq = pyfastx.Fastq(fastq_input,build_index=False)

    with open(fastq_output, 'w') as f_out:

        for name, seq, qual in fq:
            if search_func(seq):
                total_reads_found += 1
                reads_batch.append(f"@{name}\n{seq}\n+\n{qual}\n")

            if len(reads_batch) >= batch_size:
                f_out.writelines(reads_batch)
                reads_batch = []

        if reads_batch:
            f_out.writelines(reads_batch)

    logging.info(f"> Found {total_reads_found} telomere-like reads.")
    logging.info(f"> Writing telomere-like reads to: {fastq_output}")


def run_align(fastq_input,outdir,out_prefix,args):

    mapper = TelomereMapper(
        fastq_input=fastq_input,
        ref_genome=args.ref_genome,
        mm2_preset=args.mm2_preset,
        min_mapping_quality=args.min_mapq,
        out_dir=outdir,
        out_prefix=out_prefix,
        threads=args.threads
    )
    mapper.run()

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
                       default="aln",
                       help="Output prefix [default: %(default)s].")
    group.add_argument("--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")
    parser.add_argument("-t", "--threads",
                        type=int,
                        default=4,
                        help="Number of threads for minimap2 [default: %(default)s]")

    return parser

def setup_align_args(parser):

    group_name = 'Telox-align Options'
    group = parser.add_argument_group(group_name, 'Extract and align telomere-like reads')

    group.add_argument("--mm2-preset",
                       default="-ax map-ont",
                       choices=["-ax map-ont", "-ax map-pb", "-ax map-hifi"],
                       help="Minimap2 preset option [default: %(default)s]")
    group.add_argument("-q", "--min-mapq",
                       type=int,
                       default=20,
                       help="Minimum mapping quality [default: %(default)s]")

    return parser

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="telox-align: Extract and align telomere-like reads."
    )

    parser = setup_common_args(parser)
    parser = setup_align_args(parser)

    args = parser.parse_args()

    return args

def telox_align_cli():
    """
    The main entry point for the 'telox-align' command-line tool.
    """

    module = "telox-align"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    utils.setup_basic_logging()

    #### Log startup information
    utils.log_parameter_summary(args, module)

    #### Define output paths and file names
    paths = {
        'tel_reads_dir': "01_telomere_reads",
        'tel_mapping_dir': "02_telomere_mapping",
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    tel_like_fastq = os.path.join(args.tel_reads_dir,f'{args.out_prefix}.tel_like.fq')

    extract_telomeric_reads(
        fastq_input=args.long_read_fastq,
        tel_repeats=args.tel_repeats,
        fastq_output=tel_like_fastq
    )

    run_align(
        fastq_input=tel_like_fastq,
        outdir=args.outdir,
        out_prefix=args.out_prefix,
        args=args
    )

    logging.info(f"\n**** {module} completed! ****")

if __name__ == "__main__":

    telox_align_cli()

