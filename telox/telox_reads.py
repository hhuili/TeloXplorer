# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

from pathlib import Path
import gzip
import logging
import shutil
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    strand: str
    query_len: int
    query_aln_len: int
    left_softclip: int
    right_softclip: int
    mapping_category: str
    mapping_arm: str

@dataclass
class BulkConfig:
    motif: str = "TTAGGG"
    motif_C: Optional[str] = None
    motif_G: Optional[str] = None
    min_read_len: int = 1000
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
    chunk_size: int = 2000
    debug: bool = False

    @classmethod
    def from_params(cls, preset_name: Optional[str] = None, **kwargs):
        return telox_utils.apply_preset_and_overrides(cls, preset_name, PRESETS, **kwargs)

PRESETS = {
    'human': telox_utils.make_config_from_preset(BulkConfig, 'human',
        max_offset=200
    ),
    'human-r9': telox_utils.make_config_from_preset(BulkConfig, 'human-r9',
        max_offset=200
    ),
    'mouse': telox_utils.make_config_from_preset(BulkConfig, 'mouse',
        max_offset=200
    ),
    'yeast': telox_utils.make_config_from_preset(BulkConfig, 'yeast',
        max_offset=50,
        max_chimera_gap=150,
        strict_edge_repeats=3
    ),
    'arabidopsis': telox_utils.make_config_from_preset(BulkConfig, 'arabidopsis',
        max_offset=200
    ),
    'arabidopsis-r9': telox_utils.make_config_from_preset(BulkConfig, 'arabidopsis-r9',
        max_offset=200
    )
}

class BulkOutputRecord(NamedTuple):
    read_id: str
    tel_type: str
    motif_type: str
    tel_length: str
    tel_start: str
    tel_end: str
    read_length: int

def iter_fastq_reads(fastq_file, min_read_length):

    opener = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'

    with opener(fastq_file, mode) as f_fastq:
        while True:
            header = f_fastq.readline()
            if not header: break
            seq = f_fastq.readline()
            _ = f_fastq.readline() # plus
            qual = f_fastq.readline()
            if not qual:
                break

            seq = seq.strip()
            read_id = header[1:].split(' ', 1)[0].strip()

            if len(seq) >= min_read_length:
                yield read_id, seq

def classify_telomere_type(tel_len_l, tel_len_r, has_tel_l, has_tel_r):
    # mini telomere
    if has_tel_l and has_tel_r:
        return "minitel", "L,R", f"{tel_len_l},{tel_len_r}"
    # left telomere only
    elif has_tel_l:
        return "bulktel", "L", str(tel_len_l)
    # right telomere only
    elif has_tel_r:
        return "bulktel", "R", str(tel_len_r)
    # failed
    else:
        return "failed", "L,R", f"{tel_len_l},{tel_len_r}"

def select_valid_telomeres(bloom_data, aln_data, min_tel_length, min_anchor_len):

    bulktel_dict, minitel_dict, failed_dict = {}, {}, {}

    for read_id, aln_record in aln_data.items():
        read_length = aln_record.query_len

        arm_status = {}
        
        for arm in ["L", "R"]:
            tel_rec = bloom_data[arm].get(read_id)
            
            info = {
                'has_tel': False, 
                'len': 0, 
                'start': "None", 
                'end': "None",
                'raw_len': 0
            }

            if tel_rec and tel_rec.classification == "valid-TEL":

                current_subtel_len = read_length - (tel_rec.tel_len + tel_rec.initial_offset)
                is_valid_len = (tel_rec.tel_len >= min_tel_length and current_subtel_len >= min_anchor_len)

                if is_valid_len:
                    info['has_tel'] = True
                    info['len'] = tel_rec.tel_len
                    info['start'] = tel_rec.tel_start
                    info['end'] = tel_rec.tel_end

                info['raw_len'] = tel_rec.tel_len

            elif tel_rec:
                info['raw_len'] = tel_rec.classification

            arm_status[arm] = info

        has_tel_l = arm_status["L"]['has_tel']
        has_tel_r = arm_status["R"]['has_tel']
        
        len_l_group = arm_status["L"]['len'] if has_tel_l else arm_status["L"]['raw_len']
        len_r_group = arm_status["R"]['len'] if has_tel_r else arm_status["R"]['raw_len']

        telomere_type, winning_strand, winning_length = classify_telomere_type(
            len_l_group, len_r_group, has_tel_l, has_tel_r
        )

        final_start, final_end = "None", "None"
        if winning_strand == "L":
            final_start = str(arm_status["L"]['start'])
            final_end = str(arm_status["L"]['end'])
        elif winning_strand == "R":
            final_start = str(arm_status["R"]['start'])
            final_end = str(arm_status["R"]['end'])
        elif "," in winning_strand: # minitel or failed
            s_l = str(arm_status["L"]['start'])
            s_r = str(arm_status["R"]['start'])
            e_l = str(arm_status["L"]['end'])
            e_r = str(arm_status["R"]['end'])
            final_start = f"{s_l},{s_r}"
            final_end = f"{e_l},{e_r}"

        final_tel_record = BulkOutputRecord(
            read_id=read_id,
            tel_type=telomere_type,
            motif_type=winning_strand,
            tel_length=str(winning_length),
            tel_start=final_start,
            tel_end=final_end,
            read_length=read_length
        )

        target_dict = {
            'bulktel': bulktel_dict,
            'minitel': minitel_dict,
            'failed': failed_dict
        }.get(telomere_type, failed_dict)

        target_dict[read_id] = final_tel_record

    return {'bulktel': bulktel_dict,
            'minitel': minitel_dict,
            'failed': failed_dict}

def write_tel_length(tel_record_data, output):

    header = ["read_id", "tel_type", "motif_type", "tel_length", "read_length"]

    with open(output, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        for rec in tel_record_data.values():

            display_strand = rec.motif_type.replace('L', 'C-rich').replace('R', 'G-rich')

            writer.writerow([
                rec.read_id,
                rec.tel_type,
                display_strand,
                rec.tel_length,
                rec.read_length
            ])

def write_tel_seq(winning_telomeres, bloom_data, output):

    with open(output, 'w') as f_out:
        for read_id, rec in winning_telomeres.items():

            tel_rec = bloom_data.get(rec.motif_type, {}).get(read_id)
            tel_seq = tel_rec.tel_seq

            display_strand = rec.motif_type.replace('L', 'C-rich').replace('R', 'G-rich')

            header = f">{rec.read_id}|{display_strand}|{rec.tel_length}|{rec.tel_start}|{rec.tel_end}"
            f_out.write(f"{header}\n{tel_seq}\n")

def submit_chunk_tasks(executor, futures, chunk_seq_list, motif_dict, config, chunk_id):
    
    task_l = executor.submit(
        process_chunk,
        chunk_seq_list,
        motif_dict["L"],
        config.pre_merge,
        "L",
        config
    )
    futures[task_l] = "L"

    task_r = executor.submit(
        process_chunk,
        chunk_seq_list,
        motif_dict["R"],
        config.pre_merge,
        "R",
        config
    )
    futures[task_r] = "R"

def estimate_bulk_length(
        fastq_input: Union[str, Path],
        outdir: Union[str, Path],
        prefix: str,
        config: BulkConfig,
        threads: int = 4,
        logger = None
    ):

    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    motif_dict = telox_utils.parse_motifs(
        config.motif,
        getattr(config, 'motif_C', None),
        getattr(config, 'motif_G', None)
    )

    aln_data = {}
    bloom_data = {"L": {}, "R": {}}

    executor = ProcessPoolExecutor(max_workers=threads)
    futures = {}

    chunk_seq_list = []
    chunk_id = 0

    logger.info(f"Chunking reads for parallel telomere length estimation (chunk_size={config.chunk_size})")

    for read_id, seq in iter_fastq_reads(str(fastq_input), config.min_read_len):
        query_len = len(seq)

        if read_id not in aln_data:
            aln_data[read_id] = AlnRecord(
                ref_name="bulk",
                ref_start=0,
                ref_end=query_len,
                strand="+",
                query_len=query_len,
                query_aln_len=query_len,
                left_softclip=0,
                right_softclip=0,
                mapping_category="bulk",
                mapping_arm="none"
            )

            chunk_seq_list.append((read_id, "bulk", seq))

        if len(chunk_seq_list) >= config.chunk_size:
            submit_chunk_tasks(executor, futures, chunk_seq_list, motif_dict, config, chunk_id)
            chunk_seq_list = []
            chunk_id += 1

    if chunk_seq_list:
        submit_chunk_tasks(executor, futures, chunk_seq_list, motif_dict, config, chunk_id)

    total_tasks = len(futures)
    logger.info(f"Tasks submitted (n={total_tasks})")

    failed_count = 0

    for future in as_completed(futures):
        tel_arm = futures[future]
        try:
            chunk_results = future.result()
            for rec in chunk_results:
                bloom_data[tel_arm][rec.read_id] = rec
        except Exception as e:
            failed_count += 1
            logging.error(f"[ERROR] Task failed: {e}")

    executor.shutdown()

    if failed_count > 0:
        logger.info(f"WARNING: {failed_count} tasks failed. Results may be incomplete.")
    else:
        logger.info(f"All tasks completed")

    logger.info("Classifying and filtering telomere candidates")
    winning_telomeres = select_valid_telomeres(bloom_data, aln_data, config.min_tel_len, config.min_anchor_len)

    stats = {k: len(v) for k, v in winning_telomeres.items()}
    total = sum(stats.values())
    details = []

    for tel_type in ['bulktel','minitel','failed']:

        c = stats.get(tel_type, 0)
        p = (c / total * 100) if total > 0 else 0.0
        details.append(f"{tel_type}={c}({p:.1f}%)")

        output_telomere_length = outdir / f"{prefix}.{tel_type}.length.tsv"
        output_telomere_seq = outdir / f"{prefix}.{tel_type}.seq.fa"

        write_tel_length(winning_telomeres[tel_type], str(output_telomere_length))

        if tel_type in ['bulktel']:
            write_tel_seq(winning_telomeres[tel_type], bloom_data, str(output_telomere_seq))

    logger.info(f"Final telomeres: " + ", ".join(details))
