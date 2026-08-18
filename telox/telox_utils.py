# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

import logging
import subprocess
import re
import shutil
import shlex
import regex
from pathlib import Path
from dataclasses import dataclass, fields, is_dataclass
from typing import Optional, Union, List
from natsort import natsort_key

@dataclass(frozen=True)
class TeloxPreset:
    motif: Optional[str] = None
    motif_C: Optional[str] = None
    motif_G: Optional[str] = None
    baseline_density: Optional[float] = None
    fuzzy_mismatch: Optional[int] = None
    min_tel_len: Optional[int] = None
    min_anchor_len: Optional[int] = None
    bloom_opts: Optional[str] = None
    bloom_pwidth: Optional[int] = None
    bloom_ydis: Optional[float] = None
    pre_merge: Optional[str] = None
    use_pq_arms: Optional[bool] = False
    kmers: Optional[str] = None

COMMON_PRESETS = {
    'human': TeloxPreset(
        motif="TTAGGG",
        baseline_density=0.5,
        fuzzy_mismatch=2,
        min_tel_len=100,
        min_anchor_len=200,
        bloom_pwidth=200,
        bloom_ydis=0.4,
        pre_merge="no",
        use_pq_arms=True
    ),
    'human-r9': TeloxPreset(
        motif="TTAGGG|TTAAAA|CCAGGG|AAGAAG",
        motif_C="CCCTAA|CTTCTT|CCCTGG",
        motif_G="TTAGGG|TTAAAA",
        baseline_density=0.5,
        fuzzy_mismatch=2,
        min_tel_len=100,
        min_anchor_len=200,
        bloom_pwidth=200,
        bloom_ydis=0.4,
        pre_merge="no",
        use_pq_arms=True
    ),
    'mouse': TeloxPreset(
        motif="TTAGGG",
        baseline_density=0.5,
        fuzzy_mismatch=2,
        min_tel_len=100,
        min_anchor_len=200,
        bloom_pwidth=200,
        bloom_ydis=0.4,
        pre_merge="no",
        use_pq_arms=True
    ),
    'yeast': TeloxPreset(
        motif="TG{1,3}",
        baseline_density=0.7,
        fuzzy_mismatch=0,
        min_tel_len=30,
        min_anchor_len=200,
        bloom_pwidth=10,
        bloom_ydis=0.2,
        pre_merge="yes"
    ),
    'arabidopsis': TeloxPreset(
        motif="TTTAGGG",
        baseline_density=0.5,
        fuzzy_mismatch=2,
        min_tel_len=100,
        min_anchor_len=200,
        bloom_pwidth=200,
        bloom_ydis=0.4,
        pre_merge="no",
        use_pq_arms=True
    ),
    'arabidopsis-r9': TeloxPreset(
        motif="TTTAGGG|CCCAGG", 
        motif_C="CCCTAAA|CCTGGG",
        motif_G="TTTAGGG",
        baseline_density=0.5,
        fuzzy_mismatch=2,
        min_tel_len=100,
        min_anchor_len=200,
        bloom_pwidth=200,
        bloom_ydis=0.4,
        pre_merge="no",
        use_pq_arms=True
    )
}

def update_dataclass(target, source):
    if not is_dataclass(target):
        raise ValueError(f"Target {target} must be a dataclass instance")

    if is_dataclass(source):
        source_data = {f.name: getattr(source, f.name) for f in fields(source)}
    elif isinstance(source, dict):
        source_data = source
    else:
        source_data = source.__dict__ if hasattr(source, "__dict__") else {}

    target_fields = {f.name for f in fields(target)}

    for key, value in source_data.items():
        if key in target_fields and value is not None:
            setattr(target, key, value)

def make_config_from_preset(config_cls, base_name, **specifics):
    config = config_cls()
    if base_name in COMMON_PRESETS:
        update_dataclass(config, COMMON_PRESETS[base_name])
    update_dataclass(config, specifics)
    return config

def apply_preset_and_overrides(config_cls, preset_name, presets_map=None, **user_overrides):
    config = config_cls()
    target_presets = presets_map if presets_map is not None else COMMON_PRESETS
    if preset_name and target_presets:
        if preset_name not in target_presets:
            choices = ", ".join(target_presets)
            raise ValueError(
                f"Unknown preset {preset_name!r}; expected one of: {choices}."
            )
        update_dataclass(config, target_presets[preset_name])
    update_dataclass(config, user_overrides)
    return config

def simple_run_cmd(cmd, shell=False):
    if isinstance(cmd, list):
        cmd_str = shlex.join(cmd)
    else:
        cmd_str = cmd

    executable = cmd_str if shell else cmd
    
    try:
        subprocess.run(executable, shell=shell, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command failed: {cmd_str}\nStderr: {e.stderr}")

def parse_genome_indx(genome_file: Union[str, Path]):
    genome_file = Path(genome_file)
    genome_indx_file = Path(str(genome_file) + ".fai")

    if not genome_indx_file.is_file():
        logging.info(f"Genome index file not found, indexing with samtools...")
        if not genome_file.is_file():
            logging.error(f"Genome file '{genome_file}' does not exist.")
            raise FileNotFoundError(f"Genome file '{genome_file}' does not exist.")
        
        samtools_path = shutil.which('samtools')
        if not samtools_path:
            raise ValueError("Required command 'samtools' not found in PATH.")

        if genome_file.suffix == ".gz":
            cmd = f"gzip -dc {genome_file} | samtools faidx - -o {genome_indx_file}"
            subprocess.run(cmd, shell=True, check=True)
        else:
            logging.info(f"Indexing uncompressed file '{genome_file}'...")
            cmd = ['samtools', 'faidx', str(genome_file)]
            subprocess.run(cmd, check=True, capture_output=True)

    chrom_sizes = {}
    with open(genome_indx_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes

def prepare_bed_file(chrom_size_dict, output_bed, window):
    mito_regex = re.compile(r'^(chr|chromosome)?(M|MT|Mito|EBV)$', re.IGNORECASE)
    
    with open(output_bed, 'w') as f_out:
        for chrom, chrom_size in chrom_size_dict.items():
            if mito_regex.match(chrom):
                continue 

            if chrom_size <= window * 2:
                mid_point = chrom_size // 2
                if mid_point > 0:
                    f_out.write(f"{chrom}\t0\t{mid_point}\n")
                f_out.write(f"{chrom}\t{mid_point}\t{chrom_size}\n")
            else:
                f_out.write(f"{chrom}\t0\t{window}\n")
                f_out.write(f"{chrom}\t{chrom_size - window}\t{chrom_size}\n")

def reverse_complement_regex(pattern):
    complement_map = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n',
        'R': 'Y', 'Y': 'R',
        'S': 'S', 's': 's',
        'W': 'W', 'w': 'w',
        'K': 'M', 'M': 'K',
        'B': 'V', 'V': 'B',
        'D': 'H', 'H': 'D'
    }

    motifs = pattern.split('|')
    rc_motifs = []

    for motif in motifs:
        tokens = re.findall(r'(\[[^\]]+\]|\([^)]+\)|[A-Za-z])(\*|\+|\?|\{\d+,?\d*\})?', motif)
        processed_tokens = ["".join(token_tuple) for token_tuple in tokens]

        reversed_tokens = processed_tokens[::-1]

        complemented_tokens = []
        for token in reversed_tokens:
            complemented_token = "".join([complement_map.get(char, char) for char in token])
            complemented_tokens.append(complemented_token)

        rc_motifs.append("".join(complemented_tokens))

    return "|".join(rc_motifs)

def parse_motifs(motif, motif_C=None, motif_G=None):

    if motif_C and motif_G:
        return {
            "L": motif_C,
            "R": motif_G,
            "p": motif_C,
            "q": motif_G
        }

    right_arm_motif = motif
    left_arm_motif = reverse_complement_regex(right_arm_motif)
    
    return {
        "L": left_arm_motif,
        "R": right_arm_motif,
        "p": left_arm_motif,
        "q": right_arm_motif
    }

def is_dna_sequence(s):
    allowed_chars = {'A', 'T', 'C', 'G'}
    return set(s.upper()).issubset(allowed_chars)

def sequence_revcomp(s):
    s = s.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    return s[::-1]

def decode_motif_occupancy(sequence,motif):
    match_intervals = []
    occupancy = [0] * len(sequence)
    for found in re.finditer(r'{}'.format(motif),sequence):
        start = found.start(0)
        end = found.end(0)
        match_intervals.append([start,end])
        occupancy[start:end] = [1] * (end-start)

    return match_intervals, occupancy

def calc_fuzzy_density(sequence, motif, fuzzy_mismatch=2):
    fuzzy_pattern = f"({motif}){{s<={fuzzy_mismatch}}}"
    motif_found = regex.findall(fuzzy_pattern, sequence)
    found_bases = sum(len(motif) for motif in motif_found)
    seq_length = len(sequence)
    density = found_bases / seq_length if seq_length > 0 else 0.0

    return density

def roman_to_int(s):
    s = s.upper()
    roman_map = {'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100, 'D': 500, 'M': 1000}
    if not s or not all(char in roman_map for char in s):
        raise ValueError("String is not a valid Roman numeral.")
    total = 0
    prev_value = 0
    for char in reversed(s):
        value = roman_map[char]
        if value < prev_value:
            total -= value
        else:
            total += value
        prev_value = value
    if total == 0 and s != '':
         raise ValueError("Invalid Roman numeral structure.")
    return total

def get_chromosome_sort_key(chrom, organism='human'):
    if organism == 'yeast':
        try:
            return (0, roman_to_int(chrom))
        except ValueError:
            if chrom.lower() in ['m', 'mt', 'mito']:
                return (1, 0)
            else:
                return (2, natsort_key(chrom))
    else:
        return (0, natsort_key(chrom))

def create_chromosome_sorter(chrom_list):
    has_strong_roman_indicator = False
    for chrom in chrom_list:
        try:
            if chrom.isalpha() and len(chrom) > 1:
                roman_to_int(chrom)
                if chrom.lower() not in ['mt',"m","mito"]:
                    has_strong_roman_indicator = True
                    break
        except ValueError:
            pass
        if chrom.isdigit():
            has_strong_roman_indicator = False
            break
    detected_species = 'yeast' if has_strong_roman_indicator else 'human'
    return lambda c: get_chromosome_sort_key(c, organism=detected_species)

def resolve_input_samples(input_file: Optional[Path], input_list: Optional[Path]):
    samples_map = {}
    if input_list:
        input_list = Path(input_list)
        if not input_list.exists():
            raise FileNotFoundError(f"SampleSheet not found: {input_list}")
        with open(input_list, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                parts = line.split()
                if len(parts) >= 2:
                    samples_map[parts[0]] = parts[1]
                elif len(parts) == 1:
                    fpath = Path(parts[0])
                    samples_map[fpath.name.split('.')[0]] = str(fpath)
    elif input_file:
        input_file = Path(input_file)
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        fpath = Path(input_file)
        samples_map[fpath.name.split('.')[0]] = str(fpath)
    
    if not samples_map:
        raise ValueError(f"No valid samples parsed from {input_list or input_file}")

    return samples_map
