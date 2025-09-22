import os,sys
import logging
import gzip
import re
import argparse
from . import utilities as utils
from .TelomereMapper import TelomereMapper

def extract_telomeric_reads(input_fastq, tel_repeats, min_tandem_repeats, output_fastq):
    tel_repeat_dict = utils.parse_tel_repeats(tel_repeats)
    
    search_str_L = tel_repeat_dict["L"] * min_tandem_repeats
    search_str_R = tel_repeat_dict["R"] * min_tandem_repeats

    if utils.is_dna_sequence(tel_repeat_dict["L"]):
        search_func = lambda seq: search_str_L in seq or search_str_R in seq
    else:
        pattern_L = re.compile(fr"{search_str_L}")
        pattern_R = re.compile(fr"{search_str_R}")
        search_func = lambda seq: pattern_L.search(seq) or pattern_R.search(seq)

    opener = gzip.open if input_fastq.endswith('.gz') else open
    mode = 'rt' if input_fastq.endswith('.gz') else 'r'
    
    with opener(input_fastq, mode) as f_in, open(output_fastq, 'w') as f_out:
        total_reads_found = 0

        while True:
            lines = [f_in.readline() for _ in range(4)]
            if not lines[0]:
                break

            header, seq, plus_line, qual = lines
            if search_func(seq):
                f_out.write(f"{header}{seq}{plus_line}{qual}")
                total_reads_found += 1

    logging.info(f"Extracting completed. Found {total_reads_found} telomere-like reads in total.")
    logging.info(f"Results saved to {output_fastq}")


def run_telox_align(input_fastq, ref_genome, mm2_preset, min_mapping_quality, out_prefix, output_dir, threads):
    """
    telomere align.
    """
    mapper = TelomereMapper(
        out_prefix=out_prefix,
        tel_reads=input_fastq,
        ref_genome=ref_genome,
        mm2_preset=mm2_preset,
        out_dir=output_dir,
        threads=threads,
        min_mapping_quality=min_mapping_quality
    )
    mapper.run()

def parse_args():
    parser = argparse.ArgumentParser(description='Alignment Parameters')
    parser.add_argument(
        "-i", "--long-read-fastq",
        required=True,
        help="Path to the long-read FASTQ file (gzipped or plain).",
    )
    parser.add_argument(
        "-g", "--ref-genome",
        required=True,
        help="Path to the reference genome file.",
    )
    parser.add_argument(
        "-r", "--tel-repeats",
        required=True,
        help="Path to telomere repeat definition file (one repeat per line, format: arm<TAB>sequence, supports regex).",
    )
    parser.add_argument(
        "--min-tandem-repeats",
        type=int,
        default=5,
        help="Minimum number of consecutive telomeric repeat repeats (default: 5).",
    )
    parser.add_argument(
        "--minimap2-preset",
        default="-ax map-ont",
        choices=["-ax map-ont", "-ax map-pb", "-ax map-hifi"],
        help="Minimap2 preset option (default: -ax map-ont).",
    )
    parser.add_argument(
        "-q", "--min-mapping-quality",
        type=int,
        default=20,
        help="Minimum mapping quality (default: 20).",
    )
    parser.add_argument(
        "-o", "--out-prefix",
        default="telox",
        help="Output prefix (default: telox).",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: .).",
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads for mapping (default: 4).",
    )

    return parser.parse_args()

def align_cli():
    """
    The main entry point for the 'telox-align' command-line tool.
    """

    utils.setup_logging()
    args = parse_args()

    #### Define output paths and file names
    paths = {
        'tel_reads_dir': "01_telomere_reads",
        'tel_mapping_dir': "02_telomere_mapping",
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    logging.info(f"Step 1: Extracting telomere-like reads")
    tel_reads_fastq = os.path.join(args.tel_reads_dir,f'{args.out_prefix}.tel_like.fq.gz')
    extract_telomeric_reads(input_fastq=args.long_read_fastq,
                            tel_repeats=args.tel_repeats,
                            min_tandem_repeats=args.min_tandem_repeats,
                            output_fastq=tel_reads_fastq)

    logging.info(f"Step 2: Mapping telomere-like reads")
    run_telox_align(input_fastq=tel_reads_fastq,
                    ref_genome=args.ref_genome,
                    mm2_preset=args.minimap2_preset,
                    min_mapping_quality=args.min_mapping_quality,
                    out_prefix=args.out_prefix,
                    output_dir=args.tel_mapping_dir,
                    threads=args.threads)

if __name__ == "__main__":
    align_cli()

