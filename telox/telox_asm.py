# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

import logging
import csv
from pathlib import Path
from typing import Optional, Union, NamedTuple
from dataclasses import dataclass

from . import telox_utils
from . import logger as log_utils
from .telox_length import process_chunk

@dataclass(frozen=True)
class AlnRecord:
    ref_name: str
    ref_start: int
    ref_end: int
    query_seq: Optional[str]
    strand: str
    query_len: int
    query_aln_len: int
    left_softclip: int
    right_softclip: int
    mapping_category: str
    mapping_arm: str

@dataclass
class AsmConfig:
    motif: str = "TTAGGG"
    terminal_length: int = 500000
    baseline_density: float = 0.5
    pre_merge: str = "no"
    fuzzy_mismatch: int = 2
    max_offset: int = 200
    max_chimera_gap: int = 1000
    strict_edge_repeats: int = 0
    min_tel_len: int = 100
    min_anchor_len: int = 200
    bloom_pwidth: int = 200
    bloom_ydis: float = 0.4
    use_pq_arms: bool = False
    debug: bool = False

    @classmethod
    def from_params(cls, preset_name: Optional[str] = None, **kwargs):
        return telox_utils.apply_preset_and_overrides(cls, preset_name, PRESETS, **kwargs)

PRESETS = {
    'human': telox_utils.make_config_from_preset(AsmConfig, 'human', 
        terminal_length=500000,
        max_offset=-1
    ),
    'mouse': telox_utils.make_config_from_preset(AsmConfig, 'mouse', 
        terminal_length=500000,
        max_offset=-1
    ),
    'yeast': telox_utils.make_config_from_preset(AsmConfig, 'yeast', 
        terminal_length=20000,
        max_offset=100,
        max_chimera_gap=150,
        strict_edge_repeats=3
    ),
    'arabidopsis': telox_utils.make_config_from_preset(AsmConfig, 'arabidopsis', 
        terminal_length=20000,
        max_offset=-1
    )
}

class AsmOutputRecord(NamedTuple):
    query_name: str
    tel_type: str
    chrom: str
    tel_arm: str
    tel_length: str
    tel_boundary: str

def parse_assembly_file(asm_file, chrom_sizes, tolerance=10):

    aln_data = {}
    current_read_id = None
    current_chrom = None
    current_ref_start = 0
    current_ref_end = 0
    current_mapping_arm = "none"
    current_seq_parts = []

    def save_record():
        if current_read_id and current_seq_parts:
            full_seq = "".join(current_seq_parts)
            query_length = len(full_seq)
            
            aln_data[current_read_id] = AlnRecord(
                ref_name=current_chrom,
                ref_start=current_ref_start,
                ref_end=current_ref_end,
                query_seq=full_seq,
                strand="+",
                query_len=query_length,
                query_aln_len=query_length,
                left_softclip=0,
                right_softclip=0,
                mapping_category="terminal",
                mapping_arm=current_mapping_arm
            )

    with open(asm_file, 'r') as f_fasta:
        for line in f_fasta:
            line = line.strip()
            if not line: continue
            
            if line.startswith(">"):
                save_record()

                current_seq_parts = []
                current_read_id = line[1:].split(' ', 1)[0].strip()

                current_chrom = current_read_id
                current_ref_start = 0
                current_ref_end = 0
                current_mapping_arm = "none"

                if ":" in current_read_id:
                    parts = current_read_id.split(':')
                    current_chrom = parts[0]
                    try:
                        coords = parts[1].split('-')
                        current_ref_start = int(coords[0])
                        current_ref_end = int(coords[1])

                        chrom_len = chrom_sizes.get(current_chrom)
                        
                        is_left = False
                        is_right = False

                        if current_ref_start <= tolerance:
                            is_left = True

                        if chrom_len and current_ref_end >= (chrom_len - tolerance):
                            is_right = True
                        elif not chrom_len and current_ref_start > tolerance:
                            is_right = True

                        if is_left and is_right:
                            current_mapping_arm = "L,R"
                        elif is_left:
                            current_mapping_arm = "L"
                        elif is_right:
                            current_mapping_arm = "R"
                        else:
                            current_mapping_arm = "none"

                    except (IndexError, ValueError):
                        logging.warning(f"Could not parse coordinates for {current_read_id}")

                else:
                    current_mapping_arm = "L,R"

            else:
                current_seq_parts.append(line.upper())

        save_record()

    return aln_data

def classify_asm_telomere(has_tel_l, has_tel_r, mapping_arm):

    if "L" in mapping_arm and has_tel_l:
        return "passed", "L"

    if "R" in mapping_arm and has_tel_r:
        return "passed", "R"

    return "failed", "none"

def retrieve_seqs_from_aln(aln_data):
    return {read_name: aln_record.query_seq for read_name, aln_record in aln_data.items()}

def select_valid_telomeres(bloom_data, aln_data, min_tel_length, min_anchor_len):

    passed_dict, failed_dict = {}, {}

    for read_name, aln_record in aln_data.items():
        total_read_len = aln_record.query_len
        mapping_arm = aln_record.mapping_arm

        arm_status = {}
        for arm in ["L", "R"]:
            tel_rec = bloom_data[arm].get(read_name)
            info = {
                'has_tel': False,
                'len': 0, 
                'start': 0,
                'end': 0,
                'raw_len': 0
            }

            if tel_rec and tel_rec.classification == "valid-TEL":

                current_subtel_len = total_read_len - (tel_rec.tel_len + tel_rec.initial_offset)
                is_valid_len = (tel_rec.tel_len >= min_tel_length and current_subtel_len >= min_anchor_len)

                if is_valid_len:
                    info['has_tel'] = True
                    info['len'] = tel_rec.tel_len + tel_rec.initial_offset
                    info['start'] = tel_rec.tel_start
                    info['end'] = tel_rec.tel_end

                info['raw_len'] = tel_rec.tel_len

            elif tel_rec:
                info['raw_len'] = tel_rec.classification

            arm_status[arm] = info

        has_tel_l = arm_status["L"]['has_tel']
        has_tel_r = arm_status["R"]['has_tel']

        telomere_type, winning_arm = classify_asm_telomere(has_tel_l, has_tel_r, mapping_arm)

        winning_length = 0
        if winning_arm == "L":
            winning_length = arm_status["L"]['len']
        elif winning_arm == "R":
            winning_length = arm_status["R"]['len']
        else:
            # failed case: show raw lengths for debug
            winning_length = f"{arm_status['L']['raw_len']},{arm_status['R']['raw_len']}"

        tel_boundary = "NA"
        if telomere_type == "passed":
            ref_start = aln_record.ref_start
            if winning_arm == "L":
                tel_boundary = ref_start + arm_status["L"]['end'] - 1
            elif winning_arm == "R":
                tel_boundary = ref_start + arm_status["R"]['start']

        final_tel_record = AsmOutputRecord(
            query_name=read_name,
            tel_type=telomere_type,
            chrom=aln_record.ref_name,
            tel_arm=winning_arm if winning_arm != "none" else "L,R",
            tel_length=str(winning_length),
            tel_boundary=str(tel_boundary)
        )

        target_dict = {
            'passed': passed_dict,
            'failed': failed_dict
        }.get(telomere_type, failed_dict)

        target_dict[read_name] = final_tel_record

    return {'passed': passed_dict, 'failed': failed_dict}


def write_tel_length(tel_record_data, output, use_pq_arms=False):
    header = ["query_name", "tel_type", "chrom", "arm", "tel_length", "tel_boundary"]

    with open(output, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)

        for rec in tel_record_data.values():
            display_arm = rec.tel_arm.replace('L', 'p').replace('R', 'q') if use_pq_arms else rec.tel_arm
            writer.writerow([
                rec.query_name,
                rec.tel_type,
                rec.chrom,
                display_arm,
                rec.tel_length,
                rec.tel_boundary
            ])

def write_tel_seq(winning_telomeres, bloom_data, output, use_pq_arms=False):
    with open(output, 'w') as f_out:
        for read_name, rec in winning_telomeres.items():
            tel_rec = bloom_data.get(rec.tel_arm, {}).get(read_name)
            if not tel_rec:
                logging.warning(f"Sequence missing for {read_name} in {rec.tel_arm}")
                continue
            tel_seq = tel_rec.tel_seq

            display_arm = rec.tel_arm.replace('L', 'p').replace('R', 'q') if use_pq_arms else rec.tel_arm
            header = f">{rec.query_name}|{rec.chrom}|{display_arm}|{rec.tel_length}"
            f_out.write(f"{header}\n{tel_seq}\n")

def mask_assembly_fasta(input_fasta: Union[str, Path], output_fasta: Union[str, Path], passed_telomeres: dict):

    mask_coords = {}
    for rec in passed_telomeres.values():
        chrom = rec.chrom
        if chrom not in mask_coords:
            mask_coords[chrom] = {'L_end': -1, 'R_start': float('inf')}
        
        if rec.tel_boundary != "NA":
            boundary = int(rec.tel_boundary)
            if rec.tel_arm == 'L':
                mask_coords[chrom]['L_end'] = boundary
            elif rec.tel_arm == 'R':
                mask_coords[chrom]['R_start'] = boundary

    with open(input_fasta, 'r') as f_in, open(output_fasta, 'w') as f_out:
        current_chrom = None
        current_pos = 0
        L_end = -1
        R_start = float('inf')

        for line in f_in:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                f_out.write(line + "\n")

                current_chrom = line[1:].split()[0]
                current_pos = 0
                
                if current_chrom in mask_coords:
                    L_end = mask_coords[current_chrom]['L_end']
                    R_start = mask_coords[current_chrom]['R_start']
                else:
                    L_end = -1
                    R_start = float('inf')

            else:
                line_len = len(line)
                line_start = current_pos
                line_end = current_pos + line_len - 1

                if line_start > L_end and line_end < R_start:
                    f_out.write(line + "\n")

                else:
                    seq_chars = list(line)

                    if line_start <= L_end:
                        mask_stop_idx = min(line_len, L_end - line_start + 1)
                        for i in range(mask_stop_idx):
                            seq_chars[i] = 'N'

                    if line_end >= R_start:
                        mask_start_idx = max(0, R_start - line_start)
                        for i in range(mask_start_idx, line_len):
                            seq_chars[i] = 'N'

                    f_out.write("".join(seq_chars) + "\n")

                current_pos += line_len

def estimate_asm_length(
        assembly_input: Union[str, Path],
        outdir: Union[str, Path],
        prefix: str,
        config: AsmConfig,
        masked_assembly: Optional[str] = None,
        logger = None
    ):

    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    config.outdir = str(outdir)
    config.prefix = prefix
    assembly_input = Path(assembly_input)

    motif_dict = telox_utils.parse_motifs(config.motif)
    chrom_sizes = telox_utils.parse_genome_indx(str(assembly_input))

    assembly_name = assembly_input.stem
    terminal_bed_file = outdir / f'{assembly_name}.terminal.bed'
    terminal_fasta = outdir / f'{assembly_name}.terminal.fasta'

    logger.info(f"Slicing {config.terminal_length // 1000}kb sequence from chromosome ends as pseudo-long-reads")
    telox_utils.prepare_bed_file(chrom_sizes, str(terminal_bed_file), window=config.terminal_length)
    cmd = f"seqtk subseq {assembly_input} {terminal_bed_file} > {terminal_fasta}"
    logger.run_cmd(cmd, shell=True)

    bloom_data = {}

    aln_data = parse_assembly_file(str(terminal_fasta), chrom_sizes)

    for tel_arm in ["L", "R"]:

        display_arm = tel_arm.replace('L', 'p').replace('R', 'q') if config.use_pq_arms else tel_arm
        logger.info(f"Scanning {display_arm}-arm telomeres")

        motif = motif_dict[tel_arm]
        bloom_data[tel_arm] = {}

        simple_seq_list = [
            (r_name, aln_data[r_name].ref_name, aln_data[r_name].query_seq)
            for r_name in aln_data
        ]

        parsed_records = process_chunk(
            simple_seq_list=simple_seq_list,
            motif=motif,
            pre_merge=config.pre_merge,
            tel_arm=tel_arm,
            config_params=config
        )

        bloom_data[tel_arm] = {rec.read_id: rec for rec in parsed_records}

    winning_telomeres = select_valid_telomeres(bloom_data, aln_data, config.min_tel_len, config.min_anchor_len)

    stats = {k: len(v) for k, v in winning_telomeres.items()}
    total = sum(stats.values())
    details = []

    for tel_type in ['passed', 'failed']:

        c = stats.get(tel_type, 0)
        p = (c / total * 100) if total > 0 else 0.0
        details.append(f"{tel_type}={c}({p:.1f}%)")

        output_telomere_length = outdir / f"{prefix}.{tel_type}_arm.tsv"
        output_telomere_seq = outdir / f"{prefix}.{tel_type}_arm.seq.fasta"

        write_tel_length(winning_telomeres[tel_type], str(output_telomere_length), config.use_pq_arms)

        if tel_type in ['passed']:
            write_tel_seq(winning_telomeres[tel_type], bloom_data, str(output_telomere_seq), config.use_pq_arms)

    logger.info(f"Final arms: " + ", ".join(details))

    if masked_assembly and 'passed' in winning_telomeres and winning_telomeres['passed']:
        logger.info(f"Generating masked assembly -> {masked_assembly}")
        mask_assembly_fasta(assembly_input, masked_assembly, winning_telomeres['passed'])

    if not config.debug:
        if terminal_bed_file.exists():
            terminal_bed_file.unlink()
        if terminal_fasta.exists():
            terminal_fasta.unlink()
