# Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>. Licensed under GNU GPL v3.0.

import abc
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from .telox_utils import calc_fuzzy_density, decode_motif_occupancy

@dataclass
class TelomereRcord:
    read_id: str
    classification: str
    tel_start: Optional[int] = None
    tel_end: Optional[int] = None
    tel_len: Optional[int] = None
    tel_seq: Optional[str] = None
    initial_offset: Optional[int] = None
    raw_segments: List[Dict[str, Any]] = field(default_factory=list, repr=False)

class FileReader(abc.ABC):
    def __init__(self, file_path):
        self.file_path = file_path

    @abc.abstractmethod
    def __iter__(self):
        pass

    def __enter__(self):
        self._file_handle = open(self.file_path, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file_handle:
            self._file_handle.close()

class BloomParser:
    def __init__(self,
                 bloom_data,
                 seq_data,
                 motif_str,
                 baseline_density,
                 max_offset,
                 tel_arm,
                 fuzzy_mismatch,
                 max_chimera_gap=1000,
                 strict_edge_repeats=0
                 ):
        self.bloom_data = bloom_data
        self.seq_data = seq_data
        self.motif_str = motif_str
        self.baseline_density = baseline_density
        self.max_offset = max_offset
        self.tel_arm = tel_arm
        self.fuzzy_mismatch = fuzzy_mismatch
        self.max_chimera_gap = max_chimera_gap
        self.strict_edge_repeats = strict_edge_repeats
        self.match_reward = 1.0
        self.gap_open_penalty = 15.0
        self.gap_extend_penalty = 0.7

    def _refine_boundaries(self, seq: str) -> tuple[int, int]:
        if self.strict_edge_repeats <= 0:
            return 0, 0

        anchor_pattern = r'(?:' + self.motif_str + r'){' + str(self.strict_edge_repeats) + r',}'
        matches = list(re.finditer(anchor_pattern, seq.upper()))

        if not matches:
            return 0, 0

        return matches[0].start(), len(seq) - matches[-1].end()

    def _is_chimera_telomere(self, final_segments: list, seq: str) -> bool:
        seq_upper = seq.upper()
        if 'G' * 50 in seq_upper or 'C' * 50 in seq_upper:
            return True

        if self.max_chimera_gap <= 0:
            return False

        chimera_threshold = 0.3
        current_gap_len = 0
        
        for seg in final_segments:
            if seg['eff_density'] < chimera_threshold:
                current_gap_len += seg['len']
                if current_gap_len >= self.max_chimera_gap:
                    return True
            else:
                current_gap_len = 0 
                
        return False

    def _parse_and_group_segments(self):
        parsed_reads_data = defaultdict(list)
        for item in self.bloom_data:
            read_id_val = item.tchr
            parsed_reads_data[read_id_val].append({
                'start': item.start,
                'end': item.pend,
                'density': item.value,
                'id': read_id_val
            })
        return parsed_reads_data

    def __iter__(self):
        grouped_segments = self._parse_and_group_segments()

        for read_id, segments in grouped_segments.items():
            if read_id not in self.seq_data:
                yield TelomereRcord(read_id=read_id, classification="missing-sequence", raw_segments=segments)
                continue

            full_read_sequence = self.seq_data[read_id]
            read_length = len(full_read_sequence)

            sorted_segments = sorted(segments, key=lambda s: s['start'])
            labeled_segments = []

            for seg in sorted_segments:
                seg_start, seg_end = seg['start'], min(seg['end'], read_length)
                if seg_start >= seg_end:
                    continue

                seg_seq = full_read_sequence[seg_start:seg_end]
                eff_density = seg['density']

                if self.fuzzy_mismatch > 0 and eff_density < self.baseline_density:
                    fuzzy_den = calc_fuzzy_density(seg_seq, self.motif_str, self.fuzzy_mismatch)
                    eff_density = max(eff_density, fuzzy_den)

                label = 'TEL' if eff_density >= self.baseline_density else 'nonTEL'

                labeled_segments.append({
                    **seg,
                    'label': label,
                    'eff_density': eff_density,
                    'seq': seg_seq,
                    'len': seg_end - seg_start
                })

            if not labeled_segments:
                yield TelomereRcord(read_id=read_id, classification="no-valid-segments", raw_segments=segments)
                continue

            if self.tel_arm == "L":
                processing_segments = labeled_segments
            else:
                processing_segments = list(reversed(labeled_segments))

            labels = [d["label"] for d in processing_segments]
            if "TEL" not in labels:
                yield TelomereRcord(read_id=read_id, classification="no-TEL", raw_segments=labeled_segments)
                continue
            if "nonTEL" not in labels:
                yield TelomereRcord(read_id=read_id, classification="full-TEL", raw_segments=labeled_segments)
                continue

            distal_seg = processing_segments[0]
            if (self.max_offset != -1
                and 'N' not in distal_seg['seq']
                and distal_seg['label'] == "nonTEL"
                and distal_seg['len'] > self.max_offset):
                yield TelomereRcord(read_id=read_id, classification="internal-TEL", raw_segments=labeled_segments)
                continue

            start_scan_idx = -1
            for i, seg in enumerate(processing_segments):
                if seg['label'] == 'TEL':
                    start_scan_idx = i
                    break

            if start_scan_idx == -1:
                yield TelomereRcord(read_id=read_id, classification="invalid-TEL", raw_segments=labeled_segments)
                continue

            current_score = 0.0
            max_score = 0.0
            best_end_index = start_scan_idx
            baseline = self.baseline_density

            for i in range(start_scan_idx, len(processing_segments)):
                seg = processing_segments[i]
                eff_d = seg['eff_density']

                if eff_d >= baseline:
                    score_change = seg['len'] * eff_d * self.match_reward
                else:
                    is_gap_start = (i == start_scan_idx) or (processing_segments[i-1]['eff_density'] >= baseline)
                    penalty = 0.0
                    if is_gap_start:
                        penalty += self.gap_open_penalty

                    depth = baseline - eff_d
                    penalty += seg['len'] * depth * self.gap_extend_penalty
                    score_change = -penalty

                current_score += score_change

                if current_score >= max_score:
                    max_score = current_score
                    best_end_index = i

                # X-drop
                if max_score - current_score > 1000.0:
                    break

                if current_score < -50.0:
                    break

            cluster_segments = processing_segments[start_scan_idx : best_end_index + 1]
            final_segments = sorted(cluster_segments, key=lambda s: s['start'])

            tel_start = final_segments[0]['start']
            tel_end = final_segments[-1]['end']
            tel_seq = "".join(s['seq'] for s in final_segments)

            trim_l, trim_r = self._refine_boundaries(tel_seq)
            if trim_l > 0 or trim_r > 0:
                if trim_l + trim_r < len(tel_seq):
                    tel_seq = tel_seq[trim_l : len(tel_seq) - trim_r]
                    tel_start += trim_l
                    tel_end -= trim_r

            tel_span_len = tel_end - tel_start

            if self.tel_arm == "L":
                initial_offset = tel_start
            else:
                initial_offset = read_length - tel_end

            if self._is_chimera_telomere(final_segments,tel_seq):
                yield TelomereRcord(
                    read_id=read_id,
                    classification="chimera-TEL",
                    tel_start=tel_start,
                    tel_end=tel_end,
                    tel_len=tel_span_len,
                    tel_seq=tel_seq,
                    initial_offset=initial_offset,
                    raw_segments=labeled_segments
                )
                continue

            if self.max_offset != -1 and initial_offset > self.max_offset:
                yield TelomereRcord(
                    read_id=read_id, 
                    classification="internal-TEL",
                    tel_start=tel_start,
                    tel_end=tel_end,
                    tel_len=tel_span_len,
                    tel_seq=tel_seq,
                    initial_offset=initial_offset,
                    raw_segments=labeled_segments
                )
                continue

            yield TelomereRcord(
                read_id=read_id,
                classification="valid-TEL",
                tel_start=tel_start,
                tel_end=tel_end,
                tel_len=tel_span_len,
                tel_seq=tel_seq,
                initial_offset=initial_offset,
                raw_segments=labeled_segments
            )
