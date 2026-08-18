# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

from pathlib import Path
import re
import pysam
import itertools
import csv
import math
from natsort import natsorted
from typing import Optional, Union, NamedTuple
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from telox import bloom_core

from . import viz_length
from . import telox_utils
from . import logger as log_utils
from .BloomParser import BloomParser

@dataclass(frozen=True)
class AlnRecord:
    ref_name: str
    ref_start: int
    ref_end: int
    strand: str
    mapq: int
    query_len: int
    query_aln_len: int
    left_softclip: int
    right_softclip: int
    mapping_category: str
    mapping_arm: str

@dataclass
class LengthConfig:
    preset: Optional[str] = None
    motif: str = "TTAGGG"
    motif_C: Optional[str] = None
    motif_G: Optional[str] = None
    terminal_range: int = 500000
    min_read_len: int = 1000
    baseline_density: float = 0.5
    fuzzy_mismatch: int = 2
    max_offset: int = 200
    max_chimera_gap: int = 1000
    strict_edge_repeats: int = 0
    min_tel_len: int = 100
    min_anchor_len: int = 200
    bloom_pwidth: int = 200
    bloom_ydis: float = 0.4
    pre_merge: str = "no"
    chunk_size: int = 2000
    plot_length: bool = False
    debug: bool = False
    min_tel_qual: float = 0.0
    use_pq_arms: bool = False
    chrom_map: Optional[str] = None
    outdir: Optional[str] = None
    prefix: Optional[str] = None

    @classmethod
    def from_params(cls, preset_name: Optional[str] = None, **kwargs):
        return telox_utils.apply_preset_and_overrides(cls, preset_name, PRESETS, **kwargs)


def load_chromosome_map(path: Optional[Union[str, Path]]) -> dict[str, tuple[str, str]]:
    if path is None:
        return {}
    with open(path, newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        required = {"chrom", "base_chrom", "chrom_phase"}
        if not rows.fieldnames or not required.issubset(rows.fieldnames):
            raise ValueError("Chromosome map requires chrom, base_chrom, and chrom_phase columns.")
        return {row["chrom"]: (row["base_chrom"], row["chrom_phase"]) for row in rows}


def chromosome_fields(chrom: str, chrom_map: dict[str, tuple[str, str]]) -> tuple[str, str]:
    suffix_map = {
        "hap1": "hap1",
        "hap2": "hap2",
        "maternal": "mat",
        "mat": "mat",
        "paternal": "pat",
        "pat": "pat",
    }
    if chrom in chrom_map:
        base_chrom, chrom_phase = chrom_map[chrom]
        return base_chrom, suffix_map.get(str(chrom_phase).strip().lower(), chrom_phase)
    base, separator, suffix = chrom.rpartition("_")
    if separator and base and suffix.lower() in suffix_map:
        return base, suffix_map[suffix.lower()]
    return chrom, "unphased"


def display_arm(arm: str, use_pq_arms: bool = False) -> str:
    if not use_pq_arms:
        return arm
    return {"L": "p", "R": "q"}.get(arm, arm)


PRESETS = {
    'human': telox_utils.make_config_from_preset(LengthConfig, 'human',
        terminal_range=500000,
        max_offset=200,
        max_chimera_gap=1000
    ),
    'human-r9': telox_utils.make_config_from_preset(LengthConfig, 'human-r9',
        terminal_range=500000,
        max_offset=200,
        max_chimera_gap=1000
    ),
    'mouse': telox_utils.make_config_from_preset(LengthConfig, 'mouse',
        terminal_range=500000,
        max_offset=200,
        max_chimera_gap=1000
    ),
    'yeast': telox_utils.make_config_from_preset(LengthConfig, 'yeast',
        terminal_range=20000,
        max_offset=50,
        max_chimera_gap=150,
        strict_edge_repeats=3
    ),
    'arabidopsis': telox_utils.make_config_from_preset(LengthConfig, 'arabidopsis',
        terminal_range=20000,
        max_offset=200,
        max_chimera_gap=1000
    ),
    'arabidopsis-r9': telox_utils.make_config_from_preset(LengthConfig, 'arabidopsis-r9',
        terminal_range=20000,
        max_offset=200,
        max_chimera_gap=1000
    )
}

class TelomereOutputRecord(NamedTuple):
    read_id: str
    tel_type: str
    chrom: str
    tel_arm: str
    tel_length: str
    tel_start: str
    tel_end: str
    read_length: int
    mapq: int
    ref_start: str
    ref_end: str
    strand: str

def anchor_read_to_arm(ref_name, ref_start, ref_end, chrom_sizes, terminal_range, preset):

    chrom_size = chrom_sizes.get(ref_name)

    if chrom_size is None:
        return "unknown", "unknown"

    dist_to_start = ref_start
    dist_to_end = chrom_size - ref_end

    in_left_window = dist_to_start <= terminal_range
    in_right_window = dist_to_end <= terminal_range

    if in_left_window and in_right_window:
        if dist_to_start <= dist_to_end:
            return "terminal", "L"
        else:
            return "terminal", "R"
    elif in_left_window:
        return "terminal", "L"
    elif in_right_window:
        return "terminal", "R"
    else:
        return "internal", "none"

def merge_adjacent_segments(segments):
    segments.sort(key=lambda x: x[0])
    merged = []
    for seg in segments:
        if not merged or seg[0] - merged[-1][1]  > 1:
            merged.append(seg)
        else:
            merged[-1][1] = max(merged[-1][1], seg[1])

    return merged

def segment_full_sequence(tel_segments, sequence_len):
    def pairwise(iterable):
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)
    segments_flatten = [coord for seg in tel_segments for coord in seg]
    tmp = [0] + segments_flatten + [sequence_len]
    all_segments = list(pairwise(tmp))

    return all_segments

def classify_telomere_type(tel_len_l, tel_len_r, has_tel_l, has_tel_r, mapping_arm, mapping_category):

    if has_tel_l and has_tel_r:
        return "minitel", "L,R", f"{tel_len_l},{tel_len_r}"

    if has_tel_l:
        if mapping_arm == "L":
            return "chromtel", "L", str(tel_len_l)
        elif mapping_category == "internal":
            return "neotel", "L", str(tel_len_l)

    if has_tel_r:
        if mapping_arm == "R":
            return "chromtel", "R", str(tel_len_r)
        elif mapping_category == "internal":
            return "neotel", "R", str(tel_len_r)

    return "failed", "L,R", f"{tel_len_l},{tel_len_r}"

def calc_telomere_quality(qualities, start, end):
    """
    Calculates average Phred quality based on error probabilities (Log-Average).
    Formula: -10 * log10( mean(10^(-Q/10)) )
    """
    if not qualities:
        return 0.0

    q_len = len(qualities)
    start = max(0, start)
    end = min(q_len, end)

    if start >= end:
        return 0.0

    region_quals = qualities[start:end]
    count = len(region_quals)

    if count == 0:
        return 0.0

    # pysam returns integer Phred scores (already decoded from ASCII), so use q directly.
    total_error_prob = sum(math.pow(10, -q / 10.0) for q in region_quals)

    avg_error_prob = total_error_prob / count

    # convert back to Phred scale
    if avg_error_prob > 0:
        return -10.0 * math.log10(avg_error_prob)
    else:
        # theoretically impossible for standard Phred scores unless infinite quality.
        # return a safe high cap (e.g. 90).
        return 90.0

def select_valid_telomeres(bloom_data, aln_data, min_tel_length, min_anchor_len):

    chromtel_dict, neotel_dict, minitel_dict, failed_dict = {}, {}, {}, {}

    for read_id, aln_record in aln_data.items():
        chrom = aln_record.ref_name
        ref_start = aln_record.ref_start
        ref_end = aln_record.ref_end
        mapq = aln_record.mapq
        mapping_arm = aln_record.mapping_arm
        mapping_category = aln_record.mapping_category
        read_length = aln_record.query_len
        strand = aln_record.strand

        arm_status = {}
        clip_l = aln_record.left_softclip
        clip_r = aln_record.right_softclip

        for arm in ["L", "R"]:
            tel_rec = bloom_data[arm].get(read_id)
            info = {'has_tel': False, 'len': 0, 'start': "None", 'end': "None", 'raw_len': 0}

            if tel_rec and tel_rec.classification == "valid-TEL":
                distal_clip = clip_r if arm == "L" else clip_l
                current_subtel_len = read_length - (tel_rec.tel_len + tel_rec.initial_offset) - distal_clip
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
        len_l_for_grouping = arm_status["L"]['len'] if has_tel_l else arm_status["L"]['raw_len']
        len_r_for_grouping = arm_status["R"]['len'] if has_tel_r else arm_status["R"]['raw_len']

        telomere_type, winning_arm, winning_length = classify_telomere_type(
            len_l_for_grouping, len_r_for_grouping,
            has_tel_l, has_tel_r, mapping_arm, mapping_category
        )


        final_start, final_end = "None", "None"
        if winning_arm == "L":
            final_start = str(arm_status["L"]['start'])
            final_end = str(arm_status["L"]['end'])
        elif winning_arm == "R":
            final_start = str(arm_status["R"]['start'])
            final_end = str(arm_status["R"]['end'])
        elif "," in winning_arm:
            s_l, s_r = str(arm_status["L"]['start']), str(arm_status["R"]['start'])
            e_l, e_r = str(arm_status["L"]['end']), str(arm_status["R"]['end'])
            final_start = f"{s_l},{s_r}"
            final_end = f"{e_l},{e_r}"

        final_tel_record = TelomereOutputRecord(
            read_id=read_id,
            tel_type=telomere_type,
            chrom=chrom,
            tel_arm=winning_arm,
            tel_length=str(winning_length),
            tel_start=final_start,
            tel_end=final_end,
            mapq=str(mapq),
            ref_start=str(ref_start),
            ref_end=str(ref_end),
            strand=strand,
            read_length=read_length
        )

        target_dict = {
            'chromtel': chromtel_dict,
            'neotel': neotel_dict,
            'minitel': minitel_dict,
            'failed': failed_dict
        }.get(telomere_type, failed_dict)

        target_dict[read_id] = final_tel_record

    return {'chromtel': chromtel_dict,
            'neotel': neotel_dict,
            'minitel': minitel_dict,
            'failed': failed_dict}

def write_tel_length(tel_record_data, output, chrom_map=None, use_pq_arms=False):

    chrom_map = chrom_map or {}
    header = ["read_id", "tel_type", "chrom", "chrom_phase", "arm", "tel_length", "read_length", "mapq", "ref_start", "ref_end"]

    with open(output, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t', lineterminator='\n')
        writer.writerow(header)

        for rec in tel_record_data.values():
            if rec.tel_arm not in {"L", "R", "p", "q"}:
                continue
            base_chrom, chrom_phase = chromosome_fields(rec.chrom, chrom_map)
            writer.writerow([
                rec.read_id,
                rec.tel_type,
                base_chrom,
                chrom_phase,
                display_arm(rec.tel_arm, use_pq_arms),
                rec.tel_length,
                rec.read_length,
                rec.mapq,
                rec.ref_start,
                rec.ref_end
            ])

def write_tel_seq(winning_telomeres, bloom_data, qual_data, min_tel_qual, output,chrom_map=None, use_pq_arms=False):

    filtered_count = 0

    chrom_map = chrom_map or {}
    with open(output, 'w') as f_out:
        for read_id, rec in winning_telomeres.items():
            if ',' in rec.tel_start:
                continue

            final_seq = ""
            tel_rec = bloom_data.get(rec.tel_arm, {}).get(read_id)

            if tel_rec:
                if min_tel_qual > 0:
                    read_quals = qual_data.get(read_id)
                    avg_qual = calc_telomere_quality(read_quals, tel_rec.tel_start, tel_rec.tel_end)

                    if avg_qual < min_tel_qual:
                        filtered_count += 1
                        continue

                final_seq = tel_rec.tel_seq

            if final_seq:
                base_chrom, chrom_phase = chromosome_fields(rec.chrom, chrom_map)
                header = f">{rec.read_id}|{base_chrom}|{chrom_phase}|{display_arm(rec.tel_arm, use_pq_arms)}|{rec.tel_length}|{rec.tel_start}|{rec.tel_end}|{rec.strand}"
                f_out.write(f"{header}\n{final_seq}\n")

    return filtered_count

def write_tel_summary(tel_record_data, output, prefix, chrom_map=None, use_pq_arms=False):

    chrom_map = chrom_map or {}
    groups = {}
    all_lengths = []
    paternal_lengths = []
    maternal_lengths = []
    tel_type_val = None

    for record in tel_record_data.values():

        if record.tel_arm not in {"L", "R", "p", "q"}:
            continue

        if tel_type_val is None:
            tel_type_val = record.tel_type

        base_chrom, chrom_phase = chromosome_fields(record.chrom, chrom_map)
        arm = display_arm(record.tel_arm, use_pq_arms)
        length = int(record.tel_length)

        key = (base_chrom, chrom_phase, arm)

        if key not in groups:
            groups[key] = []
        groups[key].append(length)
        all_lengths.append(length)

        if chrom_phase == 'pat':
            paternal_lengths.append(length)
        elif chrom_phase == 'mat':
            maternal_lengths.append(length)

    with open(output, 'w') as f_out:
        f_out.write("level\tsample\tchrom\tchrom_phase\tarm\tcount\tTL_mean\tTL_median\n")

        if all_lengths:
            total_count = len(all_lengths)
            total_mean = sum(all_lengths) / total_count
            sorted_all = sorted(all_lengths)
            total_median = (sorted_all[total_count//2] if total_count % 2 == 1 
                          else (sorted_all[total_count//2 - 1] + sorted_all[total_count//2]) / 2)

            f_out.write(f"total\t{prefix}\tall\tall\tall\t{total_count}\t{total_mean:.0f}\t{total_median:.0f}\n")

        if paternal_lengths:
            pat_count = len(paternal_lengths)
            pat_mean = sum(paternal_lengths) / pat_count
            sorted_pat = sorted(paternal_lengths)
            pat_median = (sorted_pat[pat_count//2] if pat_count % 2 == 1 
                         else (sorted_pat[pat_count//2 - 1] + sorted_pat[pat_count//2]) / 2)

            f_out.write(f"total\t{prefix}\tpat\tpat\tall\t{pat_count}\t{pat_mean:.0f}\t{pat_median:.0f}\n")

        if maternal_lengths:
            mat_count = len(maternal_lengths)
            mat_mean = sum(maternal_lengths) / mat_count
            sorted_mat = sorted(maternal_lengths)
            mat_median = (sorted_mat[mat_count//2] if mat_count % 2 == 1 
                         else (sorted_mat[mat_count//2 - 1] + sorted_mat[mat_count//2]) / 2)

            f_out.write(f"total\t{prefix}\tmat\tmat\tall\t{mat_count}\t{mat_mean:.0f}\t{mat_median:.0f}\n")

        sorted_chroms = natsorted(groups.keys(), key=lambda x: x[0])

        for (chrom, chrom_phase, arm) in sorted_chroms:
            lengths = groups[(chrom, chrom_phase, arm)]
            count = len(lengths)
            mean_len = sum(lengths) / count if count > 0 else 0
            sorted_lengths = sorted(lengths)

            if count % 2 == 1:
                median_len = sorted_lengths[count // 2]
            else:
                median_len = (sorted_lengths[count//2 - 1] + sorted_lengths[count//2]) / 2

            f_out.write(f"chromosome\t{prefix}\t{chrom}\t{chrom_phase}\t{arm}\t{count}\t{mean_len:.0f}\t{median_len:.0f}\n")

def process_chunk(simple_seq_list, motif, pre_merge, tel_arm, config_params):

    chunk_seq_data = {item[0]: item[2] for item in simple_seq_list}
    input_data = []

    for read_id, chrom, query_seq in simple_seq_list:
        motif_segments, motif_occupancy_array = telox_utils.decode_motif_occupancy(query_seq, motif)
        motif_segments_update = segment_full_sequence(motif_segments, len(query_seq))

        if pre_merge == "yes" or sum(motif_occupancy_array) > 10000:
            motif_segments_merged = merge_adjacent_segments(motif_segments)
            motif_segments_update = segment_full_sequence(motif_segments_merged, len(query_seq))

        for segment in motif_segments_update:
            seg_start = segment[0]
            seg_end = segment[1]
            if seg_end - seg_start == 0:
                continue

            segment_density = sum(motif_occupancy_array[seg_start:seg_end]) / (seg_end - seg_start)

            input_data.append(
                bloom_core.WinDiv(chrom, read_id, seg_start, seg_end, segment_density, 1)
            )

    merged_data = []
    if input_data:
        merged_data = bloom_core.run_bloom(
            input_data=input_data,
            pwidth=config_params.bloom_pwidth,
            ydis=config_params.bloom_ydis,
            mtop=True
        )

    if config_params.debug and config_params.outdir:
        debug_file = Path(config_params.outdir) / f"{config_params.prefix}_debug_merged_{tel_arm}.txt"

        with open(debug_file, "a") as f:
            for item in merged_data:
                f.write(f"{item}\n")

    parsed_records = []
    if merged_data:
        primary_motif = motif.split('|')[0]
        bloom_parser = BloomParser(
            bloom_data=merged_data,
            seq_data=chunk_seq_data,
            motif_str=primary_motif,
            baseline_density=config_params.baseline_density,
            max_offset=config_params.max_offset,
            tel_arm=tel_arm,
            fuzzy_mismatch=config_params.fuzzy_mismatch,
            max_chimera_gap=config_params.max_chimera_gap,
            strict_edge_repeats=config_params.strict_edge_repeats
        )
        parsed_records = [record for record in bloom_parser]

    return parsed_records

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

def iter_bam_reads(bam_file: str, config: LengthConfig, chrom_sizes: dict, load_qualities: bool = False):

    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        for aln in bamfile.fetch(until_eof=True):

            if (aln.is_unmapped or
                aln.query_length < config.min_read_len or
                aln.mapping_quality < 20 or
                aln.is_secondary or
                aln.is_supplementary):
                continue

            read_id = aln.query_name
            ref_start = aln.reference_start
            ref_end = aln.reference_end

            mapping_category, mapping_arm = anchor_read_to_arm(
                aln.reference_name, ref_start, ref_end,
                chrom_sizes, config.terminal_range, config.preset
            )

            left_softclip = 0
            right_softclip = 0
            if aln.cigartuples:
                if aln.cigartuples[0][0] == 4:
                    left_softclip = aln.cigartuples[0][1]
                if aln.cigartuples[-1][0] == 4:
                    right_softclip = aln.cigartuples[-1][1]

            record = AlnRecord(
                ref_name=aln.reference_name,
                ref_start=aln.reference_start,
                ref_end=ref_end,
                strand="-" if aln.is_reverse else "+",
                mapq=aln.mapping_quality,
                query_len=aln.query_length,
                query_aln_len=aln.query_alignment_length,
                left_softclip=left_softclip,
                right_softclip=right_softclip,
                mapping_category=mapping_category,
                mapping_arm=mapping_arm
            )

            qualities = None
            if load_qualities:
                qualities = aln.query_qualities

            yield read_id, record, aln.query_sequence, qualities

def estimate_length(
        bam_file: Union[str, Path],
        ref_genome: Union[str, Path],
        outdir: Union[str, Path],
        prefix: str,
        config: LengthConfig,
        plot_config: Optional[viz_length.LengthPlotter] = None,
        threads: int = 4,
        logger=None
    ):

    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    config.outdir = str(outdir)
    config.prefix = prefix

    motif_dict = telox_utils.parse_motifs(config.motif)
    chrom_sizes = telox_utils.parse_genome_indx(ref_genome)

    aln_data = {}
    qual_data = {}
    bloom_data = {"L": {}, "R": {}}

    executor = ProcessPoolExecutor(max_workers=threads)
    futures = {}

    chunk_seq_list = []
    chunk_id = 0

    need_qual = config.min_tel_qual > 0

    logger.info(f"Chunking reads for parallel telomere length estimation (chunk_size={config.chunk_size})")

    for read_id, record, sequence, qualities in iter_bam_reads(str(bam_file), config, chrom_sizes, load_qualities=need_qual):

        if read_id not in aln_data:
            aln_data[read_id] = record
            if need_qual:
                qual_data[read_id] = qualities
            chunk_seq_list.append((read_id, record.ref_name, sequence))

        if len(chunk_seq_list) >= config.chunk_size:
            submit_chunk_tasks(executor, futures, chunk_seq_list, motif_dict, config, chunk_id)
            chunk_seq_list = []
            chunk_id += 1

    if chunk_seq_list:
        submit_chunk_tasks(executor, futures, chunk_seq_list, motif_dict, config, chunk_id)
        chunk_id += 1

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
            logger.info(f"[ERROR] Task failed (Arm: {tel_arm}): {e}")

    executor.shutdown()

    if failed_count > 0:
        logger.info(f"WARNING: {failed_count} tasks failed. Results may be incomplete.")
    else:
        logger.info(f"All tasks completed")

    logger.info(f"Classifying and filtering telomere candidates")
    chrom_map = load_chromosome_map(config.chrom_map)

    winning_telomeres = select_valid_telomeres(
        bloom_data,
        aln_data,
        config.min_tel_len,
        config.min_anchor_len
    )

    stats = {k: len(v) for k, v in winning_telomeres.items()}
    total = sum(stats.values())
    details = []

    for tel_type in ['chromtel','neotel','minitel','failed']:

        c = stats.get(tel_type, 0)
        p = (c / total * 100) if total > 0 else 0.0
        details.append(f"{tel_type}={c}({p:.1f}%)")

        output_telomere_length = outdir / f"{prefix}.{tel_type}.length.tsv"
        output_telomere_seq = outdir / f"{prefix}.{tel_type}.seq.fasta"

        write_tel_length(
            winning_telomeres[tel_type],
            str(output_telomere_length),
            chrom_map,
            use_pq_arms=config.use_pq_arms,
        )

        if tel_type in ['chromtel','neotel']:

            filtered_count = write_tel_seq(
                winning_telomeres[tel_type],
                bloom_data,
                qual_data,
                config.min_tel_qual,
                str(output_telomere_seq),
                chrom_map,
                use_pq_arms=config.use_pq_arms,
            )

    logger.info(f"Final telomeres: " + ", ".join(details))

    for tel_type in ['chromtel']:
        output_telomere_summary = outdir / f"{prefix}.{tel_type}.summary.tsv"
        write_tel_summary(
            winning_telomeres[tel_type],
            str(output_telomere_summary),
            prefix,
            chrom_map,
            use_pq_arms=config.use_pq_arms,
        )

        if config.plot_length:
            length_file = outdir / f"{prefix}.{tel_type}.length.tsv"
            if plot_config is None:
                plot_config = viz_length.LengthPlotter()

            logger.info(f"Plotting chromosome-end telomere length distributions ({tel_type})")
            viz_length.plot_length(
                length_input=str(length_file),
                outdir=outdir,
                prefix=f'{prefix}.{tel_type}',
                config=plot_config
            )
