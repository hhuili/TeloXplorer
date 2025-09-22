import subprocess
import logging
import shutil
import re
import regex
import os
import sys

class CustomFormatter(logging.Formatter):
    def format(self, record):
        if getattr(record, 'notime', False):
            return record.getMessage()
        else:
            return super().format(record)

def setup_logging():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    # Apply custom formatter to the root handler
    handler = logging.getLogger().handlers[0]
    custom_fmt = CustomFormatter(handler.formatter._fmt, datefmt=handler.formatter.datefmt)
    handler.setFormatter(custom_fmt)

def validate_command_exists(command_name):
    if not shutil.which(command_name):
        raise ValueError(f"Required command '{command_name}' not found in PATH.")

def run_cmd(command, use_shell=False):

    validate_command_exists(command[0])

    cmd_str = ' '.join(command) if isinstance(command, (list, tuple)) else command

    logging.info(f"Running the command: {cmd_str}")
    try:
        subprocess.run(cmd_str if use_shell else command, shell=use_shell, check=True)
    except subprocess.CalledProcessError as e:
        logging.info(f"Error: {e}")
        raise

def parse_sample_table(input_file):
    sample_table = {}
    with open(input_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.strip().split('\t')
            sample_table[line[0]] = line[1]

    return sample_table

def parse_chromsize(genome_indx_file):

    genome_file = genome_indx_file.replace(".fai","")
    if not os.path.isfile(genome_indx_file):
        logging.info("no genome index file find, run samtools faidx...")
        cmd = ['samtools','faidx',genome_file]
        run_cmd(cmd)

    chromsize_dict = {}
    with open(genome_indx_file) as f:
        for line in f:
            line = line.strip().split('\t')
            chromsize_dict[line[0]] = int(line[1])

    return chromsize_dict

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

def prepare_chrom_window_file(chromsize_dict,output_window,window_size=500000,middle_window="no"):
    with open(output_window,'w') as f_out:
        for chrm,chrm_size in chromsize_dict.items():
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




