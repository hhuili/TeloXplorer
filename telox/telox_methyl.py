# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

from __future__ import annotations

import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Optional, Set, Tuple

import pysam

from . import logger as log_utils


MM_HEADER_RE = re.compile(r"^([ACGTUN])([+-])([a-z]+|[0-9]+)([.?])?$")
CHEBI_ALIASES = {
    "27551": "m",  # 5mC
    "76792": "h",  # 5hmC
}


@dataclass
class MethylConfig:
    c_threshold: float = 0.66
    mod_m_threshold: float = 0.90
    mod_h_threshold: float = 0.90
    chunk_size: int = 1000
    subtel_extension: int = 10000
    valid_mod_codes: Set[str] = field(default_factory=lambda: {"m", "h"})
    canonical_base: str = "C"
    check_canonical: bool = True
    debug: bool = False

    @classmethod
    def from_params(cls, **kwargs):
        config = cls()
        for key, value in kwargs.items():
            if hasattr(config, key) and value is not None:
                setattr(config, key, value)
        config.canonical_base = str(config.canonical_base).upper()
        config.valid_mod_codes = {str(code) for code in config.valid_mod_codes}
        return config


def relative_telomere_position(eff_pos: int, tel_start: int, tel_length: int, arm: str,) -> int:
    p_local = eff_pos - tel_start

    if arm in {"L", "p"}:
        return tel_length - 1 - p_local
    if arm in {"R", "q"}:
        return p_local
    raise ValueError(f"Unsupported telomere arm: {arm!r}")


def extract_read_info(fasta_path: Path, id_output_path: Path):
    read_info = {}
    with open(fasta_path) as f_in, open(id_output_path, "w") as f_ids:
        for line in f_in:
            if not line.startswith(">"):
                continue

            parts = line.strip()[1:].split("|")
            if len(parts) < 8:
                continue

            strand = parts[-1]
            tel_start = int(parts[-3])
            tel_length = int(parts[-4])
            tel_arm = parts[-5]
            chrom_phase = parts[-6]
            chrom = parts[-7]
            read_id = "|".join(parts[:-7])

            read_info[read_id] = {
                "chrom": chrom,
                "chrom_phase": chrom_phase,
                "arm": tel_arm,
                "tel_len": tel_length,
                "tel_start": tel_start,
                "strand": strand,
            }
            f_ids.write(f"{read_id}\n")

    return read_info


def filter_modbam(
    modbam: Path,
    id_list: Path,
    output_bam: Path,
    threads: int = 4,
    logger=None,
):
    if output_bam.exists():
        logger.info(f"BAM exists, overwriting: {output_bam}")

    logger.run_cmd(
        [
            "samtools",
            "view",
            "-b",
            "-N",
            str(id_list),
            "-@",
            str(threads),
            "-o",
            str(output_bam),
            str(modbam),
        ]
    )


def normalize_mod_code(code: object) -> str:
    """Normalize a pysam/SAM modification code to the Telox short code."""
    text = str(code)
    return CHEBI_ALIASES.get(text, text)


def parse_mm_skip_modes(mm_tag: str) -> Dict[str, Dict[str, str]]:
    """Return ``base -> modification -> skip mode`` from an MM tag.

    ``.`` (and an omitted flag) means skipped bases may be treated as having
    low modification probability. ``?`` means the skipped bases have unknown
    modification status. Multiple MM groups for the same canonical base are
    kept separately by modification code instead of overwriting one another.
    """
    modes: DefaultDict[str, Dict[str, str]] = defaultdict(dict)

    for block in mm_tag.split(";"):
        if not block:
            continue
        header = block.split(",", 1)[0]
        match = MM_HEADER_RE.fullmatch(header)
        if match is None:
            continue

        base, _strand, raw_codes, mode = match.groups()
        mode = mode or "."
        codes = list(raw_codes) if raw_codes.isalpha() else [raw_codes]

        for raw_code in codes:
            code = normalize_mod_code(raw_code)
            previous = modes[base].get(code)
            # If duplicate declarations disagree, '?' is the conservative
            # interpretation: an omitted probability cannot be assumed zero.
            modes[base][code] = "?" if "?" in {previous, mode} else "."

    return dict(modes)


def ml_probability(quality: int) -> Optional[float]:
    if quality is None or not 0 <= quality <= 255:
        return None
    return (quality + 0.5) / 256.0


def classify_call(
    probs: Dict[str, float],
    canonical_base: str,
    canonical_threshold: float,
    mod_thresholds: Dict[str, float],
) -> Optional[Tuple[str, float]]:
    p_canonical = max(0.0, 1.0 - sum(probs.values()))
    best_state = canonical_base
    best_prob = p_canonical

    for mod_code, probability in probs.items():
        if probability > best_prob:
            best_state = mod_code
            best_prob = probability

    threshold = (
        canonical_threshold
        if best_state == canonical_base
        else mod_thresholds.get(best_state, canonical_threshold)
    )
    return (best_state, best_prob) if best_prob >= threshold else None


def parse_and_filter_mods(
    read_data: Dict[str, Any],
    meta: Dict[str, Any],
    config: MethylConfig,
):
    mod_dict = read_data.get("modified_bases") or {}
    query_sequence = read_data.get("query_sequence") or ""
    seq_len = read_data.get("query_length", 0)
    mm_tag = read_data.get("MM_tag") or ""

    if seq_len == 0 or not query_sequence or not mm_tag:
        return []

    tel_start = meta["tel_start"]
    tel_len = meta["tel_len"]
    arm = meta["arm"]
    strand = meta["strand"]
    ext_sub = config.subtel_extension

    target_eff_start = 0
    target_eff_end = seq_len - 1
    if arm in {"L", "p"}:
        target_eff_end = tel_start + tel_len + ext_sub - 1
    elif arm in {"R", "q"}:
        target_eff_start = tel_start - ext_sub
    else:
        raise ValueError(f"Unsupported telomere arm: {arm!r}")

    if strand == "+":
        start_idx = max(0, target_eff_start)
        end_idx = min(seq_len, target_eff_end + 1)
    elif strand == "-":
        start_idx = max(0, seq_len - 1 - target_eff_end)
        end_idx = min(seq_len, seq_len - target_eff_start)
    else:
        raise ValueError(f"Unsupported read strand: {strand!r}")

    if start_idx >= end_idx:
        return []

    mm_skip_modes = parse_mm_skip_modes(mm_tag)
    canonical_base = config.canonical_base
    base_modes = mm_skip_modes.get(canonical_base)
    if not base_modes:
        return []

    site_probs: DefaultDict[int, Dict[str, float]] = defaultdict(dict)
    for base_tuple, sites in mod_dict.items():
        if len(base_tuple) < 3:
            continue

        base_name = str(base_tuple[0]).upper()
        mod_code = normalize_mod_code(base_tuple[2])
        if base_name != canonical_base:
            continue

        for pos, quality in sites:
            if not start_idx <= pos < end_idx:
                continue
            if config.check_canonical and query_sequence[pos].upper() != canonical_base:
                continue

            probability = ml_probability(quality)
            if probability is not None:
                site_probs[pos][mod_code] = probability

    mod_thresholds = {
        "m": config.mod_m_threshold,
        "h": config.mod_h_threshold,
    }
    valid_mods = []

    for pos in range(start_idx, end_idx):
        if query_sequence[pos].upper() != canonical_base:
            continue

        probs = site_probs.get(pos, {})

        # For a '?' MM group, absence of an explicit probability means the
        # modification state is unknown, so the full C/m/h state vector cannot
        # be reconstructed safely at this position.
        if any(mode == "?" and code not in probs for code, mode in base_modes.items()):
            continue

        result = classify_call(
            probs,
            canonical_base=canonical_base,
            canonical_threshold=config.c_threshold,
            mod_thresholds=mod_thresholds,
        )
        if result is None:
            continue

        state, probability = result
        if state != canonical_base and state not in config.valid_mod_codes:
            continue

        eff_pos = pos if strand == "+" else seq_len - 1 - pos
        p_rel = relative_telomere_position(eff_pos, tel_start, tel_len, arm)
        valid_mods.append((p_rel, state, probability))

    valid_mods.sort(key=lambda item: item[0])
    return valid_mods


def worker_process_batch(
    tasks: List[Tuple[Dict[str, Any], Dict[str, Any], MethylConfig]],
):
    out_lines = []

    for read_data, meta, config in tasks:
        valid_mods = parse_and_filter_mods(read_data, meta, config)
        if not valid_mods:
            continue

        m_sites, h_sites, valid_sites = [], [], []
        for p_rel, mod_type, _probability in valid_mods:
            site = str(int(p_rel))
            valid_sites.append(site)
            if mod_type == "m":
                m_sites.append(site)
            elif mod_type == "h":
                h_sites.append(site)

        out_lines.append(
            "\t".join(
                [
                    read_data["query_name"],
                    meta["chrom"],
                    meta["chrom_phase"],
                    meta["arm"],
                    str(meta["tel_len"]),
                    str(read_data["query_length"]),
                    ",".join(m_sites) if m_sites else ".",
                    ",".join(h_sites) if h_sites else ".",
                    ",".join(valid_sites),
                ]
            )
            + "\n"
        )

    return out_lines


def find_methyl(
    modbam: Path,
    telseq: Path,
    outdir: Path,
    prefix: str,
    config: MethylConfig,
    threads: int = 4,
    logger=None,
):
    logger = logger or log_utils.get_logger()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    output_id_list = outdir / f"{prefix}.read_ids.txt"
    tel_bam = outdir / f"{prefix}.tel.bam"
    output_plot = outdir / f"{prefix}.mods.plot_data.tsv"

    read_info = extract_read_info(telseq, output_id_list)
    if not read_info:
        logger.warning(f"No valid telomere headers found in {telseq}. Exiting!")
        return

    logger.info("Extracting target telomere reads from modBAM")
    filter_modbam(modbam, output_id_list, tel_bam, threads=threads, logger=logger)

    save_verbosity = pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(str(tel_bam), "rb", check_sq=False)
    pysam.set_verbosity(save_verbosity)

    task_buffer = []
    logger.info("Parsing modifications and mapping relative to T-S boundaries")

    with open(output_plot, "w") as f_out, ProcessPoolExecutor(max_workers=threads) as executor:
        f_out.write(
            "read_id\tchrom\tchrom_phase\tarm\ttel_length\tread_length\t"
            "mods_m\tmods_h\tvalid_sites\n"
        )

        futures = []
        target_ids = set(read_info)

        for read in bam:
            if read.query_name not in target_ids:
                continue
            if read.is_secondary or read.is_supplementary:
                continue

            query_sequence = read.query_sequence or ""
            read_len = len(query_sequence) if query_sequence else (read.query_length or 0)
            if read_len == 0:
                continue

            try:
                mm_tag = read.get_tag("MM")
            except KeyError:
                try:
                    mm_tag = read.get_tag("Mm")
                except KeyError:
                    continue

            # SAM recommends MN as a consistency check after operations such as
            # hard clipping. If present and inconsistent, MM/ML coordinates are
            # stale and must not be interpreted.
            try:
                mn_tag = int(read.get_tag("MN"))
            except KeyError:
                mn_tag = None
            if mn_tag is not None and mn_tag != read_len:
                logger.warning(
                    f"Skipping {read.query_name}: MN={mn_tag} but SEQ length={read_len}."
                )
                continue

            modified_bases = read.modified_bases or {}
            raw_read_data = {
                "query_name": read.query_name,
                "query_sequence": query_sequence,
                "query_length": read_len,
                "modified_bases": dict(modified_bases),
                "MM_tag": mm_tag,
            }
            task_buffer.append((raw_read_data, read_info[read.query_name], config))

            if len(task_buffer) >= config.chunk_size:
                futures.append(executor.submit(worker_process_batch, task_buffer))
                task_buffer = []

                done = []
                for index, future in enumerate(futures):
                    if not future.done():
                        continue
                    try:
                        f_out.writelines(future.result())
                    except Exception as exc:
                        logger.error(f"Worker error: {exc}")
                    done.append(index)
                for index in reversed(done):
                    del futures[index]

        if task_buffer:
            futures.append(executor.submit(worker_process_batch, task_buffer))

        for future in as_completed(futures):
            try:
                f_out.writelines(future.result())
            except Exception as exc:
                logger.error(f"Final flush error: {exc}")

    bam.close()

    if not config.debug:
        for temporary_path in (tel_bam, output_id_list):
            if temporary_path.exists():
                temporary_path.unlink()
