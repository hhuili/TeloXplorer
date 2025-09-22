import os,sys
import logging
import shutil
import pysam
import argparse
import itertools
from pathlib import Path
from . import utilities as utils
from .BloomParser import BloomParser

PRESETS = {
    'Human': {
        'min_telomere_freq': 0.6,
        'max_repeat_mismatches': 2,
        'min_initial_telomere_offset': 10000,
        'min_telomere_length': 100,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': False,
    },
    'Mouse': {
        'min_telomere_freq': 0.6,
        'max_repeat_mismatches': 2,
        'min_initial_telomere_offset': 10000,
        'min_telomere_length': 100,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': False,
    },
    'Yeast': {
        'min_telomere_freq': 0.7,
        'max_repeat_mismatches': 0,
        'min_initial_telomere_offset': 10000,
        'min_telomere_length': 30,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 50 --ydis 0.2",
        'primary_merge': True,
    }
}

def add_mapping_arm(ref_name, ref_start, ref_end, chromosome_sizes, subtelomere_boundary):
    if ref_name not in chromosome_sizes:
        chrom_group = "unknown_chrom"
        mapping_arm = "unknown_chrom"
        return chrom_group, mapping_arm

    chrom_size = chromosome_sizes[ref_name]
    if ref_start <= subtelomere_boundary:
        chrom_group = "chromTerminal"
        mapping_arm = "L"
    elif ref_end >= chrom_size - subtelomere_boundary:
        chrom_group = "chromTerminal"
        mapping_arm = "R"
    else:
        chrom_group = "chromInternal"
        mapping_arm = "none"
    return chrom_group, mapping_arm

def parse_assembly_file(asm_file):

    aln_data = {}
    with open(asm_file,'r') as f_fasta:
        for line in f_fasta:
            if line.startswith(">"):
                read_id = line[1:].split(' ', 1)[0].strip()
                if ":" in read_id:
                    chrom = read_id.split(':')[0]
                    ref_start = int(read_id.split(':')[1].split('-')[0])
                    ref_end = int(read_id.split(':')[1].split('-')[1])
                    mapping_arm = "L" if ref_start == 1 else "R"
                else:
                    chrom = read_id
                    ref_start = 0
            else:
                query_seq = line.strip().upper()
                query_length = len(query_seq)

                aln_data[read_id] = {'ref_name': chrom,
                                     'ref_start': ref_start,
                                     'ref_end': ref_end,
                                     'query_seq': query_seq,
                                     'strand': "+",
                                     'query_len':query_length,
                                     'query_align_len': query_length,
                                     'left_softclip': 0,
                                     'right_softclip': 0,
                                     'chrom_group': "chromTerminal",
                                     'mapping_arm': mapping_arm}

        return aln_data

def parse_alignment_file(aln_file,min_mapq,min_read_length,chrom_sizes,chrom_subtel_boundary):

    with pysam.AlignmentFile(aln_file, "rb") as bamfile:
        aln_data = {}
        for aln in bamfile.fetch(until_eof=True):
            if (aln.is_unmapped or
                    aln.mapping_quality < min_mapq or
                    aln.query_length < min_read_length or
                    aln.is_secondary or
                    aln.is_supplementary):
                continue

            read_id = aln.query_name
            ref_start = aln.reference_start
            ref_end = aln.reference_start+aln.reference_length
            chrom_group, mapping_arm = add_mapping_arm(aln.reference_name,
                                                       ref_start,
                                                       ref_end,
                                                       chrom_sizes,
                                                       chrom_subtel_boundary)

            # Extract soft-clipped lengths from CIGAR tuples
            left_softclip = 0
            right_softclip = 0
            if aln.cigartuples:
                # CIGAR operation for soft-clip is 4
                if aln.cigartuples[0][0] == 4:
                    left_softclip = aln.cigartuples[0][1]
                if aln.cigartuples[-1][0] == 4:
                    right_softclip = aln.cigartuples[-1][1]

            if read_id not in aln_data:
                aln_data[read_id] = {'ref_name':aln.reference_name,
                                        'ref_start':aln.reference_start,
                                        'ref_end':ref_end,
                                        'query_seq':aln.query_sequence,
                                        'strand':"-" if aln.is_reverse else "+",
                                        'query_len':aln.query_length,
                                        'query_align_len':aln.query_alignment_length,
                                        'left_softclip':left_softclip,
                                        'right_softclip':right_softclip,
                                        'chrom_group':chrom_group,
                                        'mapping_arm':mapping_arm}
            else:
                print(read_id)

    return aln_data

def load_sequences_from_aln_reader(aln_data):
    return {read_name: aln_record['query_seq'] for read_name,aln_record in aln_data.items()}

def merge_adjacent_tel_segments(segments):
    segments.sort(key=lambda x: x[0])

    merged = []
    for seg in segments:
        # If the list is empty, or the current segment does not coincide with the previous segment, add it directly
        if not merged or seg[0] - merged[-1][1]  > 1:
            merged.append(seg)
        else:
            # Otherwise, we can merge with the previous seg
            merged[-1][1] = max(merged[-1][1], seg[1])

    return merged

def add_nontel_segments(segments,sequence):

    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = itertools.tee(iterable)
        next(b, None)
        return zip(a, b)

    ## add nontel segments in tel_segments
    segments_flatten = sum(segments,[])
    tmp = [0,len(sequence)]
    tmp[1:-1] = segments_flatten

    segments_update = list(pairwise(tmp))

    return segments_update

def prepare_tel_segments(aln_data,tel_repeat,output_file,primary_merge=False):

    with open(output_file,'w') as f_out:
        for read_name in aln_data:
            chrom = aln_data[read_name]['ref_name']
            query_seq = aln_data[read_name]['query_seq']
            tel_segments, tel_repeat_freq = utils.calculate_exact_telomere_freq(query_seq,tel_repeat)
            tel_segments_update = add_nontel_segments(tel_segments,query_seq)
            if primary_merge or sum(tel_repeat_freq) > 10000:
                tel_segments_merged = merge_adjacent_tel_segments(tel_segments)
                tel_segments_update = add_nontel_segments(tel_segments_merged,query_seq)
            for segment in tel_segments_update:
                seg_start = segment[0]
                seg_end = segment[1]
                if seg_end - seg_start == 0:
                    continue
                mean_tel_repeat_freq = sum(tel_repeat_freq[seg_start:seg_end])/(seg_end-seg_start)
                f_out.write(f"{chrom}\t{read_name}\t{seg_start}\t{seg_end}\t{mean_tel_repeat_freq}\n")


def assign_tel_group(tel_length_l,tel_length_r,has_tel_l,has_tel_r,mapping_arm,chrom_group):

    if has_tel_l and has_tel_r:
        tel_group = f'{chrom_group}_miniTel'
        winning_tel_arm = "L,R"
        winning_tel_length = f'{tel_length_l},{tel_length_r}'
    elif has_tel_l and (mapping_arm == "L" or chrom_group == "chromInternal"):
        tel_group = chrom_group
        winning_tel_arm = "L"
        winning_tel_length = tel_length_l
    elif has_tel_r and (mapping_arm == "R" or chrom_group == "chromInternal"):
        tel_group = chrom_group
        winning_tel_arm = "R"
        winning_tel_length = tel_length_r
    else:
        tel_group = 'filterTel'
        winning_tel_arm = "L,R"
        winning_tel_length = f'{tel_length_l},{tel_length_r}'

    return tel_group,winning_tel_arm,winning_tel_length

def get_winning_telomeres(bloom_data,aln_data,min_tel_length,min_subtel_length):

    chromtel_dict,neotel_dict,minitel_dict,filtertel_dict = {},{},{},{}
    for read_name,aln_record in aln_data.items():
        tel_record_l = bloom_data["L"][read_name]
        tel_record_r = bloom_data["R"][read_name]
        chrom = aln_record['ref_name']
        mapping_start = aln_record['ref_start']
        query_align_length = aln_record['query_align_len']
        left_softclip = aln_record['left_softclip']
        right_softclip = aln_record['right_softclip']
        mapping_end = aln_record['ref_end']
        mapping_arm = aln_record['mapping_arm']
        chrom_group = aln_record['chrom_group']

        if tel_record_l.classification == "TEL-positive":
            tel_length_l = tel_record_l.tel_len + tel_record_l.initial_tel_offset
            subtel_length_l = query_align_length + left_softclip - tel_length_l
            has_tel_l = True if tel_length_l >= min_tel_length and subtel_length_l >= min_subtel_length else False
        else:
            tel_length_l = tel_record_l.classification
            subtel_length_l = query_align_length
            has_tel_l = False
        if tel_record_r.classification == "TEL-positive":
            tel_length_r = tel_record_r.tel_len + tel_record_r.initial_tel_offset
            subtel_length_r = query_align_length + right_softclip - tel_length_r
            has_tel_r = True if tel_length_r >= min_tel_length and  subtel_length_r >= min_subtel_length else False
        else:
            tel_length_r = tel_record_r.classification
            subtel_length_r = query_align_length
            has_tel_r = False

        tel_group,winning_tel_arm,winning_tel_length = assign_tel_group(tel_length_l,tel_length_r,has_tel_l,has_tel_r,mapping_arm,chrom_group)

        if tel_group == "chromTerminal":
            tel_boundary = winning_tel_length if winning_tel_arm == "L" else int(mapping_end) - int(winning_tel_length)
        else:
            tel_boundary = "NA"

        final_tel_record = [read_name,
                            tel_group,
                            chrom,
                            winning_tel_arm,
                            str(winning_tel_length),
                            str(tel_boundary)]

        if tel_group == "filterTel":
            # filtered telomeres
            filtertel_dict[read_name] = final_tel_record
        elif tel_group == "chromTerminal":
            ## chromosome-end telomeres
            chromtel_dict[read_name] = final_tel_record
        elif tel_group == "chromInternal":
            ## neo-telomeres (chromosome-internal telomeres)
            neotel_dict[read_name] = final_tel_record
        else:
            ## mini telomeres (both end telomeres)
            minitel_dict[read_name] = final_tel_record

    winning_telomeres = {'chromtel': chromtel_dict,
                       'neotel': neotel_dict,
                       'minitel': minitel_dict,
                       'filtertel': filtertel_dict}

    return winning_telomeres

def run_bloom(tel_segment_file,bloom_options,output_dir,bloom_output,tel_arm):

    bloom_outdir = os.path.join(output_dir,f"bloom_{tel_arm}")
    if not os.path.isdir(bloom_outdir):
        os.makedirs(bloom_outdir)

    bloom_cmd = ["BLOOM",bloom_options,"--mtop","--input",tel_segment_file,"--outdir",bloom_outdir]
    utils.run_cmd(bloom_cmd,use_shell=True)

    cmd = ["cp",f'{bloom_outdir}/merge.final.*.out',bloom_output]
    utils.run_cmd(cmd,use_shell=True)

def write_tel_length(tel_record_data,output):

    with open(output,'w') as f_out:
        f_out.write(f'Query_name\tTel_group\tChrom\tTel_arm\tTel_length\tTel_boundary\n')
        for tel_record in tel_record_data.values():
            f_out.write('\t'.join(tel_record)+'\n')

def write_tel_seq(winning_telomeres,bloom_data,output):

    with open(output,'w') as f_out:
        for read_name, tel_record in winning_telomeres.items():
            chrom = tel_record[2]
            tel_arm = tel_record[3]
            tel_length = tel_record[4]
            tel_seq = bloom_data[tel_arm][read_name].tel_seq
            f_out.write(f'>{read_name}|{chrom}|{tel_arm}|{tel_length}\n{tel_seq}\n')

def run_asm_length(asm_file, tel_repeats, bloom_options, out_prefix, outdir, args):
    """
    Calculating assembly telomere length.
    """

    #### arguments
    min_telomere_freq=args.min_telomere_freq
    min_initial_telomere_offset=args.min_initial_telomere_offset
    max_repeat_mismatches=args.max_repeat_mismatches
    min_telomere_length=args.min_telomere_length
    min_subtelomere_length=args.min_subtelomere_length
    primary_merge=args.primary_merge

    #### create chromosome size dict
    tel_repeat_dict = utils.parse_tel_repeats(tel_repeats)
    chromsize_dict = utils.parse_chromsize(f"{asm_file}.fai")

    #### prepare terminal 500k fasta file of the assembly
    path_obj = Path(asm_file)
    asm_file_name = path_obj.stem
    terminal_500k_window = os.path.join(outdir,f'{asm_file_name}.500k_window.bed')
    genome_terminal_500k_fasta = os.path.join(outdir,f'{asm_file_name}.500k.fa')

    utils.prepare_chrom_window_file(chromsize_dict,terminal_500k_window,500000)
    cmd = ['seqtk','subseq',asm_file,terminal_500k_window,'>',genome_terminal_500k_fasta]
    utils.run_cmd(cmd,use_shell=True)

    bloom_data = {}
    for tel_arm in ["L", "R"]:
        tel_repeat = tel_repeat_dict[tel_arm]
        bloom_data[tel_arm] = {}

        # Paths for intermediate files
        out_tel_segments = os.path.join(outdir, f'{out_prefix}.tel_segments_{tel_arm}.txt')
        bloom_output = os.path.join(outdir, f'{out_prefix}.bloom_{tel_arm}.txt')

        aln_data = parse_assembly_file(genome_terminal_500k_fasta)
        prepare_tel_segments(aln_data, tel_repeat, out_tel_segments, primary_merge=primary_merge)
        run_bloom(out_tel_segments, bloom_options, outdir, bloom_output, tel_arm)

        seq_data = load_sequences_from_aln_reader(aln_data)

        bloom_parser = BloomParser(
            bloom_file=bloom_output,
            seq_data=seq_data,
            tel_repeat_str=tel_repeat,
            min_telomere_labeling_threshold=min_telomere_freq,
            min_initial_telomere_offset=min_initial_telomere_offset,
            tel_arm=tel_arm,
            fuzzy_freq_threshold=min_telomere_freq,
            max_repeat_mismatches=max_repeat_mismatches
        )
        tel_records = {tel_record.read_id: tel_record for tel_record in bloom_parser}
        bloom_data[tel_arm] = tel_records

        ## remove intermediate file
        if not args.debug:
            os.remove(out_tel_segments)
            os.remove(bloom_output)
            shutil.rmtree(f'{outdir}/bloom_{tel_arm}')

    ## calculate telomere length and filtering
    winning_telomeres = get_winning_telomeres(bloom_data,aln_data,min_telomere_length,min_subtelomere_length)

    ## write final telomere reads
    tel_groups = ['chromtel','neotel','minitel','filtertel']
    logging.info(f"output telomere length...\n")
    for tel_group in tel_groups:
        output_telomere_length = os.path.join(outdir,f"{out_prefix}.{tel_group}_length.tsv")
        output_telomere_seq = os.path.join(outdir,f"{out_prefix}.{tel_group}.tel_seq.fa")
        write_tel_length(winning_telomeres[tel_group],output_telomere_length)
        if tel_group in ['chromtel','neotel']:
            write_tel_seq(winning_telomeres[tel_group],bloom_data,output_telomere_seq)

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="telox_asm: Teloxplorer module for assembly telomere length analysis."
    )
    parser.add_argument(
        "--preset",
        type=str,
        choices=['Human', 'Mouse', 'Yeast'],
        required=False,
        help="Use a preset configuration for a specific species (e.g., Human, Mouse, Yeast). "
            "Manually specified parameters will override preset values."
    )
    parser.add_argument(
        "-i", "--asm-fasta",
        required=True,
        help="Input assembly file (fasta).",
    )
    parser.add_argument(
        "-r", "--tel-repeats",
        required=True,
        help="Path to telomere repeat definition file (one repeat per line, format: arm<TAB>sequence, supports regex).",
    )
    parser.add_argument(
        "--min-telomere-freq",
        type=float,
        default=0.6,
        help="Initial scan: frequency threshold for distinguishing telomere from non-telomere segments (default: 0.6).",
    )
    parser.add_argument(
        "--primary-merge",
        action='store_true',
        required=False,
        help="Merge adjacent telomere segments before running BLOOM. (Default: False).",
    )
    parser.add_argument(
        "--max-repeat-mismatches",
        type=int,
        default=2,
        help="Secondary scan: mismatch tolerance for re-labeling non-telomere segments within telomere segments (default: 2).",
    )
    parser.add_argument(
        "--min-initial-telomere-offset",
        type=int,
        default=200,
        help="Minimum non-telomere sequence at the head of a telomere read (default: 200 bp).",
    )
    parser.add_argument(
        "--min-telomere-length",
        type=int,
        default=100,
        help="Minimum telomere length for a telomere read (default: 100 bp).",
    )
    parser.add_argument(
        "--min-subtelomere-length",
        type=int,
        default=100,
        help="Minimum sub-telomere length for a telomere read (default: 100 bp).",
    )
    parser.add_argument(
        "--bloom-options",
        type=str,
        default="--bloom_pwidth 200 --bloom_ydis 0.4",
        help="BLOOM preset options (default: --bloom_pwidth 200 --bloom_ydis 0.4)."
            "Used for telomere segments merge.",
    )
    parser.add_argument(
        "-o", "--out-prefix",
        default="telox",
        help="Output prefix.",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory.",
    )
    parser.add_argument(
        "--debug",
        action='store_true',
        required=False,
        help="Enable debug mode for verbose logging.",
    )

    args = parser.parse_args()
    user_provided_args = set()
    for action in parser._actions:
        if action.dest != 'help' and hasattr(args, action.dest):
            if getattr(args, action.dest) != action.default:
                user_provided_args.add(action.dest)

    args.user_provided_args = user_provided_args

    return args


def asm_cli():
    """
    The main entry point for the 'telox-asm' command-line tool.
    """

    utils.setup_logging()
    args = parse_arguments()

    if args.preset:
        logging.info(f"Applying preset for '{args.preset}':")
        preset_values = PRESETS.get(args.preset)
        if not preset_values:
            available_presets = list(PRESETS.keys())
            raise ValueError(f"Invalid preset '{args.preset}'. Available presets are: {available_presets}")

        for key, value in preset_values.items():
            if key not in args.user_provided_args and hasattr(args, key):
                setattr(args, key, value)
                logging.info(f"  -> Setting '{key}' to preset value: {value}", extra={'notime': True})

    #### Define output paths and file names
    paths = {
        'tel_length_dir': "telomere_length"
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    logging.info(f"Calculating assembly telomere length")
    run_asm_length(asm_file=args.asm_fasta,
                    tel_repeats=args.tel_repeats,
                    bloom_options=args.bloom_options,
                    out_prefix=args.out_prefix,
                    outdir=args.tel_length_dir,
                    args=args)

if __name__ == "__main__":
    asm_cli()

