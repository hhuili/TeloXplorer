import abc
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from .utilities import calculate_fuzzy_telomere_freq


@dataclass
class TelomereRcord:
    """
    A structured data container for the telomere analysis result of a single read.
    """
    read_id: str
    classification: str # e.g., 'non-TEL', 'all-TEL', 'internal-TEL', 'TEL-positive'
    
    # Optional fields, only populated for TEL-positive reads
    tel_start: Optional[int] = None
    tel_end: Optional[int] = None
    tel_len: Optional[int] = None
    tel_seq: Optional[str] = None
    initial_tel_offset: Optional[int] = None
    
    # Internal data for debugging or advanced analysis
    raw_segments: List[Dict[str, Any]] = field(default_factory=list, repr=False)

class FileReader(abc.ABC):
    """
    Abstract Base Class for file parsers that read genomic data.
    Subclasses should implement the 'parse' method as a generator.
    """
    def __init__(self, file_path):
        self.file_path = file_path

    @abc.abstractmethod
    def __iter__(self):
        """Yields records from the file."""
        pass

    def __enter__(self):
        self._file_handle = open(self.file_path, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file_handle:
            self._file_handle.close()

class BloomParser(FileReader):
    """
    Parses output from a bloom filter, identifies telomeric regions in reads,
    and calculates telomere lengths.

    This version processes the entire bloom file upon initialization and allows
    for direct, on-demand access to the analysis result for each read ID.
    """
    def __init__(self,
                 bloom_file,
                 seq_data,
                 tel_repeat_str,
                 min_labeling_threshold,
                 max_initial_offset,
                 tel_arm,
                 fuzzy_freq_threshold,
                 max_mismatch):
        """
        Initializes the BloomParser.

        Args:
            bloom_file (str): Path to the bloom filter output file.
            seq_data (dict): Input sequences for each read (dict).
            tel_repeat_str (str): The telomere repeat string (e.g., 'TTAGGG').
            min_telomere_labeling_threshold (float): Frequency threshold to label a segment as 'TEL'.
            max_initial_telomere_offset (int): Minimum length for a non-terminal 'nonTEL' segment.
            tel_arm (str): The telomere arm to process ('L' or 'R').
            fuzzy_freq_threshold (float): Frequency threshold for fuzzy matching.
            max_repeat_mismatch (int): Max allowed mismatches for fuzzy repeat search.
        """
        super().__init__(bloom_file)
        
        self.seq_data = seq_data
        self.tel_repeat_str = tel_repeat_str
        self.min_labeling_threshold = min_labeling_threshold
        self.max_initial_offset = max_initial_offset
        self.tel_arm = tel_arm
        self.fuzzy_freq_threshold = fuzzy_freq_threshold
        self.max_mismatch = max_mismatch

    def _parse_and_group_segments(self):
        """Reads the bloom file and groups segments by read ID."""
        parsed_reads_data = defaultdict(list)
        with open(self.file_path, 'r') as f_in:
            for line_num, line_content in enumerate(f_in, 1):
                parts = line_content.strip().split('\t')
                if len(parts) < 5:
                    continue
                
                read_id_val = parts[1]
                try:
                    segment_start_val = int(parts[2])
                    segment_end_val = int(parts[3])
                    tel_freq_val = float(parts[4])
                except ValueError:
                    logging.warning(f"Skipping line {line_num} due to invalid numeric value: {line_content.strip()}")
                    continue
                
                parsed_reads_data[read_id_val].append({
                    'start': segment_start_val, 'end': segment_end_val,
                    'tel_freq': tel_freq_val, 'id': read_id_val
                })

        return parsed_reads_data

    def _merge_adjacent_tel_segments(self, segments_list):
        """Merges adjacent 'TEL' labeled segments into a single segment."""
        if not segments_list:
            return []
        merged_output = []
        idx = 0
        while idx < len(segments_list):
            segment = segments_list[idx]
            if segment['label'] == 'TEL':
                tel_block = [segment]
                next_idx = idx + 1
                while next_idx < len(segments_list) and segments_list[next_idx]['label'] == 'TEL':
                    tel_block.append(segments_list[next_idx])
                    next_idx += 1
                
                block_start = min(s['start'] for s in tel_block)
                block_end = max(s['end'] for s in tel_block)
                sorted_segments = sorted(tel_block, key=lambda s: s['start'])
                block_sequence = "".join(s['seq'] for s in sorted_segments)

                merged_output.append({
                    'id': segment['id'], 'start': block_start, 'end': block_end,
                    'tel_freq': -1, 'label': 'TEL', 'len': block_end - block_start,
                    'seq': block_sequence
                })
                idx = next_idx
            else:
                merged_output.append(segment)
                idx += 1
        return merged_output

    def __iter__(self):

        grouped_segments = self._parse_and_group_segments()

        for read_id, segments in grouped_segments.items():

            if read_id not in self.seq_data:
                yield TelomereRcord(read_id=read_id, classification="sequence-missing", raw_segments=segments)
                continue

            full_read_sequence = self.seq_data[read_id]
            read_length = len(full_read_sequence)
            
            # 1. Initial Labeling & Canonical Ordering
            sorted_segments = sorted(segments, key=lambda s: s['start'])
            labeled_segments = []
            for seg in sorted_segments:
                seg_start, seg_end = seg['start'], min(seg['end'], read_length)
                if seg_start >= seg_end: continue
                label = 'TEL' if seg['tel_freq'] >= self.min_labeling_threshold else 'nonTEL'
                labeled_segments.append({
                    **seg, 'label': label,
                    'seq': full_read_sequence[seg_start:seg_end],
                    'len': seg_end - seg_start
                })

            if not labeled_segments:
                yield TelomereRcord(read_id=read_id, classification="no-valid-segments", raw_segments=segments)
                continue

            # 2. Arm-specific processing order
            processing_segments = labeled_segments if self.tel_arm == "L" else list(reversed(labeled_segments))

            # 3. Case Analysis & Classification
            labels = [d["label"] for d in processing_segments]
            if "TEL" not in labels:
                yield TelomereRcord(read_id=read_id, classification="non-TEL", raw_segments=labeled_segments)
                continue
            if "nonTEL" not in labels:
                yield TelomereRcord(read_id=read_id, classification="all-TEL", raw_segments=labeled_segments)
                continue
            
            num_segments = len(processing_segments)
            if num_segments == 2 and processing_segments[-1]['label'] == 'TEL':
                yield TelomereRcord(read_id=read_id, classification="all-TEL", raw_segments=labeled_segments)
                continue
            if (self.max_initial_offset != -1 
                and 'N' not in processing_segments[0]['seq']
                and processing_segments[0]['label'] == "nonTEL" 
                and processing_segments[0]['len'] > self.max_initial_offset):
                yield TelomereRcord(read_id=read_id, classification="internal-TEL", raw_segments=labeled_segments)
                continue

            # Perfect Cases
            if num_segments in [2, 3] and processing_segments[1]['label'] == 'TEL':
                initial_tel_offset = processing_segments[0]['len']
                res_seg = processing_segments[1]
                yield TelomereRcord(
                    read_id=read_id,
                    classification="TEL-positive",
                    tel_start=res_seg['start'],
                    tel_end=res_seg['end'],
                    tel_len=res_seg['len'],
                    tel_seq=res_seg['seq'],
                    initial_tel_offset=initial_tel_offset,
                    raw_segments=labeled_segments
                )
                continue

            # 4. Iterative Relabeling for Irregular Patterns
            working_segments = [dict(s) for s in processing_segments]
            if self.max_mismatch > 0:
                for _ in range(3): # Max iterations
                    changed = False
                    indices_to_relabel = [i for i, s in enumerate(working_segments) if 0 < i < len(working_segments) - 1 and s['label'] == 'nonTEL']
                    if not indices_to_relabel:
                        break

                    for i in indices_to_relabel:
                        tel_freq_fuzzy = calculate_fuzzy_telomere_freq(working_segments[i]['seq'], self.tel_repeat_str, self.max_mismatch)
                        if tel_freq_fuzzy >= self.fuzzy_freq_threshold:
                            working_segments[i]['label'] = 'TEL'
                            changed = True
                    
                    if changed:
                        working_segments = self._merge_adjacent_tel_segments(working_segments)
                    else:
                        break
            
            # 5. Final Telomere Length Calculation
            found_tel = False
            for i, seg_data in enumerate(working_segments):
                if i <= 1 and seg_data['label'] == 'TEL':
                    if i == 0:
                        initial_tel_offset = 0
                    else:
                        initial_tel_offset = working_segments[0]['len']

                    yield TelomereRcord(
                        read_id=read_id, classification="TEL-positive", tel_start=seg_data['start'], tel_end=seg_data['end'],
                        tel_len=seg_data['len'], tel_seq=seg_data['seq'], initial_tel_offset=initial_tel_offset,
                        raw_segments=labeled_segments
                    )
                    found_tel = True
                    break

            if not found_tel:
                yield TelomereRcord(read_id=read_id, classification="TEL-negative", raw_segments=labeled_segments)



