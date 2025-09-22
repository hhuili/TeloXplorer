import os,sys
import io
import logging
import shutil
import pysam
import argparse
import itertools
import subprocess
from . import utilities as utils
from .BloomParser import BloomParser

PRESETS = {
    'Human': {
        'min_tandem_repeats': 5,
        'chrom_subtelomere_boundary': 500000,
        'min_read_length': 1000,
        'min_mapping_quality': 20,
        'min_telomere_freq': 0.6,
        'max_repeat_mismatches': 2,
        'min_initial_telomere_offset': 200,
        'min_telomere_length': 100,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': False,
    },
    'Mouse': {
        'min_tandem_repeats': 5,
        'chrom_subtelomere_boundary': 500000,
        'min_read_length': 1000,
        'min_mapping_quality': 20,
        'min_telomere_freq': 0.6,
        'max_repeat_mismatches': 2,
        'min_initial_telomere_offset': 200,
        'min_telomere_length': 100,
        'min_subtelomere_length': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': False,
    },
    'Yeast': {
        'min_tandem_repeats': 5,
        'chrom_subtelomere_boundary': 20000,
        'min_read_length': 1000,
        'min_mapping_quality': 20,
        'min_telomere_freq': 0.7,
        'max_repeat_mismatches': 0,
        'min_initial_telomere_offset': 100,
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
        return

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
                #logging.info(f'Duplicated read: {read_id}')
                pass

    return aln_data

def load_sequences_from_aln_reader(aln_data):
    return {read_name: aln_record['query_seq'] for read_name,aln_record in aln_data.items()}

def merge_adjacent_tel_segments(segments):
    segments.sort(key=lambda x: x[0])

    merged = []
    for seg in segments:
        # If the list is empty, or the current segment does not coincide with the previous segment, add it directly
        if not merged or seg[0] - merged[-1][1]  > 10:
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

def prepare_tel_segments_generator(aln_data, tel_repeat):

    for read_name in aln_data:
        chrom = aln_data[read_name]['ref_name']
        query_seq = aln_data[read_name]['query_seq']
        tel_segments, tel_repeat_freq = utils.calculate_exact_telomere_freq(query_seq, tel_repeat)
        tel_segments_update = add_nontel_segments(tel_segments, query_seq)
        for segment in tel_segments_update:
            seg_start = segment[0]
            seg_end = segment[1]
            if seg_end - seg_start == 0:
                continue
            mean_tel_repeat_freq = sum(tel_repeat_freq[seg_start:seg_end]) / (seg_end - seg_start)

            yield f"{chrom}\t{read_name}\t{seg_start}\t{seg_end}\t{mean_tel_repeat_freq}\n"

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

def check_telomere_homopolymer(telomere_seq,threshold=50):

    sequence_upper = telomere_seq.upper()
    g_pattern = 'G' * (threshold + 1)
    c_pattern = 'C' * (threshold + 1)

    return g_pattern in sequence_upper or c_pattern in sequence_upper

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
            tel_seq_l = tel_record_l.tel_seq
            tel_length_l = tel_record_l.tel_len
            subtel_length_l = query_align_length + left_softclip - (tel_length_l + tel_record_l.initial_tel_offset)
            has_tel_l = (tel_length_l >= min_tel_length and subtel_length_l >= min_subtel_length)
            if has_tel_l:
                has_tel_l = not check_telomere_homopolymer(tel_seq_l)
        else:
            tel_length_l = tel_record_l.classification
            subtel_length_l = query_align_length
            has_tel_l = False
        if tel_record_r.classification == "TEL-positive":
            tel_seq_r = tel_record_r.tel_seq
            tel_length_r = tel_record_r.tel_len
            subtel_length_r = query_align_length + right_softclip - (tel_length_r + tel_record_r.initial_tel_offset)
            has_tel_r = (tel_length_r >= min_tel_length and subtel_length_r >= min_subtel_length)
            if has_tel_r:
                has_tel_r = not check_telomere_homopolymer(tel_seq_r)
        else:
            tel_length_r = tel_record_r.classification
            subtel_length_r = query_align_length
            has_tel_r = False

        tel_group,winning_tel_arm,winning_tel_length = assign_tel_group(tel_length_l,tel_length_r,has_tel_l,has_tel_r,mapping_arm,chrom_group)

        final_tel_record = [read_name,
                            tel_group,
                            chrom,
                            winning_tel_arm,
                            str(winning_tel_length),
                            str(mapping_start),
                            str(mapping_end)]

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

def run_bloom_in_memory(segment_generator, bloom_options):
    """
    Run BLOOM in memory.
    - segment_generator: A generator that yields input data lines.
    - bloom_options: Command line arguments for BLOOM.
    
    Returns: BLOOM's standard output content (in string format).
    """

    bloom_cmd = ["BLOOM", bloom_options, "--mtop", "--input", "-"]

    # Use Popen to start the process and connect its stdin, stdout, and stderr streams
    process = subprocess.Popen(
        bloom_cmd,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=False
    )

    # Write data from the generator line by line to BLOOM's stdin
    for line in segment_generator:
        process.stdin.write(line)

    # Close stdin to indicate that data transmission is complete
    process.stdin.close()
    
    # Read all output from stdout and stderr
    stdout_data, stderr_data = process.communicate()
    
    # Check if any error occurred
    if process.returncode != 0:
        print("BLOOM process failed with error:")
        print(stderr_data)
        return None

    return stdout_data

def write_tel_length(tel_record_data,output):

    with open(output,'w') as f_out:
        f_out.write(f'Query_name\tTel_group\tChrom\tTel_arm\tTel_length\tMapping_start\tMapping_end\n')
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

def run_telox_length(aln_file, tel_repeats, ref_genome, bloom_options, out_prefix, outdir, args):
    """
    Calculating telomere length.
    """

    #### arguments
    min_mapq=args.min_mapping_quality
    min_read_length=args.min_read_length
    min_telomere_freq=args.min_telomere_freq
    min_initial_telomere_offset=args.min_initial_telomere_offset
    chrom_subtelomere_boundary=args.chrom_subtelomere_boundary
    max_repeat_mismatches=args.max_repeat_mismatches
    min_telomere_length=args.min_telomere_length
    min_subtelomere_length=args.min_subtelomere_length
    primary_merge=args.primary_merge

    tel_repeat_dict = utils.parse_tel_repeats(tel_repeats)
    chromsize_dict = utils.parse_chromsize(f"{ref_genome}.fai")

    bloom_data = {}
    for tel_arm in ["L", "R"]:
        tel_repeat = tel_repeat_dict[tel_arm]
        bloom_data[tel_arm] = {}

        # Paths for intermediate files
        out_tel_segments = os.path.join('/dev/shm', f'{out_prefix}.tel_segments_{tel_arm}.txt')
        bloom_output = os.path.join(outdir, f'{out_prefix}.bloom_{tel_arm}.txt')

        aln_data = parse_alignment_file(aln_file,min_mapq,min_read_length,chromsize_dict,chrom_subtelomere_boundary)
        prepare_tel_segments(aln_data, tel_repeat, out_tel_segments, primary_merge=primary_merge)
        run_bloom(out_tel_segments, bloom_options, outdir, bloom_output, tel_arm)
        #segment_generator = prepare_tel_segments_generator(aln_data, tel_repeat)
        #bloom_output_str = run_bloom_in_memory(segment_generator, bloom_options)
        seq_data = load_sequences_from_aln_reader(aln_data)
        #bloom_output_stream = io.StringIO(bloom_output_str)
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
        if os.path.exists(out_tel_segments):
            os.remove(out_tel_segments)
        if not args.debug:
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
        description="telox_length: Teloxplorer module for chromosome-specific telomere length analysis."
    )
    parser.add_argument(
        "--preset",
        type=str,
        choices=['Human', 'Mouse', 'Yeast'],
        required=False,
        help="Use a preset configuration for a specific species (e.g., Human, Mouse, Yeast). "
            "Manually specified parameters will override preset values."
    )

    # Group 1: Input Files
    input_group = parser.add_argument_group("Input Files", "Input files for the pipeline.")
    input_group.add_argument(
        "-i", "--aln-file",
        required=True,
        help="Input alignment file (bam/sam).",
    )
    input_group.add_argument(
        "-g", "--ref-genome",
        required=True,
        help="Path to the reference genome file.",
    )
    input_group.add_argument(
        "-r", "--tel-repeats",
        required=True,
        help="Path to telomere repeat definition file (one repeat per line, format: arm<TAB>sequence, supports regex).",
    )

    # Group 2: Telomere Length Parameters
    telomere_group = parser.add_argument_group("Telomere Length Parameters", "Parameters for calculating telomere length.")
    telomere_group.add_argument(
        "-q", "--min-mapping-quality",
        type=int,
        default=20,
        help="Minimum mapping quality (default: 20).",
    )
    telomere_group.add_argument(
        "--chrom-subtelomere-boundary",
        type=int,
        default=500000,
        help="Telomere-subtelomere boundary of chromosomes (default: 500000 bp).",
    )
    telomere_group.add_argument(
        "--min-read-length",
        type=int,
        default=1000,
        help="Minimum read length for telomere analysis (default: 1000 bp).",
    )
    telomere_group.add_argument(
        "--min-telomere-freq",
        type=float,
        default=0.6,
        help="Initial scan: frequency threshold for distinguishing telomere from non-telomere segments (default: 0.6).",
    )
    telomere_group.add_argument(
        "--primary-merge",
        action='store_true',
        required=False,
        help="Merge adjacent telomere segments before running BLOOM. (Default: False).",
    )
    telomere_group.add_argument(
        "--max-repeat-mismatches",
        type=int,
        default=2,
        help="Secondary scan: mismatch tolerance for re-labeling non-telomere segments within telomere segments (default: 2).",
    )
    telomere_group.add_argument(
        "--min-initial-telomere-offset",
        type=int,
        default=200,
        help="Minimum non-telomere sequence at the head of a telomere read (default: 200 bp).",
    )
    telomere_group.add_argument(
        "--min-telomere-length",
        type=int,
        default=100,
        help="Minimum telomere length for a telomere read (default: 100 bp).",
    )
    telomere_group.add_argument(
        "--min-subtelomere-length",
        type=int,
        default=100,
        help="Minimum sub-telomere length for a telomere read (default: 100 bp).",
    )

    # Group 3: BLOOM Parameters
    bloom_group = parser.add_argument_group("BLOOM Parameters", "Options for the BLOOM tool.")
    bloom_group.add_argument(
        "--bloom-options",
        type=str,
        default="--bloom_pwidth 200 --bloom_ydis 0.4",
        help="BLOOM preset options (default: --bloom_pwidth 200 --bloom_ydis 0.4)."
            "Used for telomere segments merge.",
    )

    # Group 4: Output and Debugging
    output_group = parser.add_argument_group("Output and Debugging", "Output directory and debug settings.")
    output_group.add_argument(
        "-o", "--out-prefix",
        default="telox",
        help="Output prefix (default: telox).",
    )
    output_group.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: .).",
    )
    output_group.add_argument(
        "--debug",
        action="store_true",
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


def length_cli():
    """
    The main entry point for the 'telox-length' command-line tool.
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

    logging.info(f"Calculating telomere length")
    run_telox_length(aln_file=args.aln_file,
                    tel_repeats=args.tel_repeats,
                    ref_genome=args.ref_genome,
                    bloom_options=args.bloom_options,
                    out_prefix=args.out_prefix,
                    outdir=args.tel_length_dir,
                    args=args)

if __name__ == "__main__":
    length_cli()

