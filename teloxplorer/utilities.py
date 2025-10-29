import subprocess
import logging
import shutil
import re
import regex
import os
import sys

def setup_logging(log_file, debug=False):

    log_level = logging.DEBUG if debug else logging.INFO

    logger = logging.getLogger()
    logger.setLevel(log_level)

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    formatter = logging.Formatter('%(message)s')

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

def setup_basic_logging(debug=False):

    log_level = logging.DEBUG if debug else logging.INFO

    logging.basicConfig(
        level=log_level,
        format='%(message)s',
        stream=sys.stdout,
        force=True
    )

def log_parameter_summary(args, module):

    logging.info(f"{module} parameters: ")

    if hasattr(args,"preset"):
        logging.info(f"  --preset: {args.preset}")

    if hasattr(args, "long_read_fastq"):
        logging.info(f"  --long-read-fastq: {args.long_read_fastq}")
    if hasattr(args, "assembly"):
        logging.info(f"  --assembly: {args.assembly}")
    if hasattr(args, "tel_fasta"):
        logging.info(f"  --tel-fasta: {args.tel_fasta}")
    if hasattr(args, "bed_methyl"):
        logging.info(f"  --bed-methyl: {args.bed_methyl}")
    if hasattr(args, "ref_genome"):
        logging.info(f"  --ref-genome: {args.ref_genome}")
    if hasattr(args, "genome_subtel_range"):
        logging.info(f"  --genome-subtel-range: {args.genome_subtel_range}")
    if hasattr(args, "tel_repeats"):
        tel_repeat_dict = parse_tel_repeats(args.tel_repeats)
        logging.info(f"  --tel-repeats: {tel_repeat_dict}")

    if hasattr(args, "min_tel_freq"):
        logging.info(f"  --min-tel-freq: {args.min_tel_freq}")
    if hasattr(args, "mm2_preset"):
        logging.info(f"  --mm2-preset: {args.mm2_preset}")
    if hasattr(args, "min_mapq"):
        logging.info(f"  --min-mapq: {args.min_mapq}")
    if hasattr(args, "bloom_options"):
        logging.info(f"  --bloom-options: {args.bloom_options}")

    if hasattr(args, "out_prefix"):
        logging.info(f"  --out-prefix: {args.out_prefix}")
    if hasattr(args, "outdir"):
        logging.info(f"  --outdir: {args.outdir}")
    if hasattr(args, "threads"):
        logging.info(f"  --threads: {args.threads}")
    
    if hasattr(args, "find_TVR"):
        logging.info(f"  --find-TVR: {args.find_TVR}")
        if args.find_TVR == "yes":
            logging.info(f"  --kmers: {args.kmers}")

def validate_command_exists(command_name):

    if not shutil.which(command_name):
        raise ValueError(f"Required command '{command_name}' not found in PATH.")

def run_cmd(command, use_shell=False):

    validate_command_exists(command[0])
    cmd_str = ' '.join(command) if isinstance(command, (list, tuple)) else command

    logging.info(f"> Running the command: {cmd_str}")
    subprocess.run(cmd_str if use_shell else command, shell=use_shell, check=True)

def parse_genome_indx(genome_indx_file):

    genome_file = genome_indx_file.replace(".fai","")
    if not os.path.isfile(genome_indx_file):
        logging.info("> Genome index not found, generating with samtools faidx...")
        cmd = ['samtools','faidx',genome_file]
        run_cmd(cmd)

    chrom_sizes = {}
    with open(genome_indx_file) as f:
        for line in f:
            line = line.strip().split('\t')
            chrom_sizes[line[0]] = int(line[1])

    return chrom_sizes

def parse_tel_repeats(repeat_file):

    tel_repeat_dict = {}
    with open(repeat_file) as f:
        for line in f:
            line = line.strip().split('\t')
            tel_repeat_dict[line[0]] = line[1]

    return tel_repeat_dict

def add_path(args, name, path):

    if not os.path.isdir(path):
        os.makedirs(path)
    setattr(args, name, path)
    return args

def natural_sort_key(chrom):
    """natural sorting of chromosome names"""

    chrom_clean = chrom.replace('chr', '').split('_')[0]
    special_order = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
    
    if chrom_clean.isdigit():
        return (int(chrom_clean), chrom)
    elif chrom_clean in special_order:
        return (special_order[chrom_clean], chrom)
    return (999, chrom)


def prepare_window_file(chrom_sizes,output_window,window_size=500000,middle_window="no"):

    with open(output_window,'w') as f_out:
        for chrm,chrm_size in chrom_sizes.items():
            if chrm.lower() not in ["chrm", "chrmito", "chrebv"]:
                ## start window
                window_start,window_end = 0,window_size
                f_out.write(f"{chrm}\t{window_start}\t{window_end}\n")
                ## middle window
                if middle_window == "yes":
                    window_start = window_size
                    window_end = chrm_size - window_size
                    f_out.write(f"{chrm}\t{window_start}\t{window_end}\n")
                ## end window
                window_start,window_end = chrm_size - window_size,chrm_size
                f_out.write(f"{chrm}\t{window_start}\t{window_end}\n")
            else:
                window_start = '0'
                f_out.write(f"{chrm}\t{window_start}\t{chrm_size}\n")

def is_dna_sequence(s):

    allowed_chars = {'A', 'T', 'C', 'G'}
    return set(s.upper()).issubset(allowed_chars)

def sequence_revcomp(s):

    s = s.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    return s[::-1]

def calculate_exact_telomere_freq(sequence,repeat_pattern):

    tel_segments = []
    tel_repeat_freq = [0] * len(sequence)

    for repeat_find in re.finditer(r'{}'.format(repeat_pattern),sequence):
        start = repeat_find.start(0)
        end = repeat_find.end(0)
        tel_segments.append([start,end])
        tel_repeat_freq[start:end] = [1] * (end-start)

    return tel_segments, tel_repeat_freq


def calculate_fuzzy_telomere_freq(sequence, tel_repeat, max_mismatches=2):

    fuzzy_pattern = f"({tel_repeat}){{s<={max_mismatches}}}"

    repeat_found = regex.findall(fuzzy_pattern, sequence)
    repeat_bases = sum(len(repeat) for repeat in repeat_found)

    seq_len = len(sequence)
    frequency = repeat_bases / seq_len if seq_len > 0 else 0.0

    return frequency



