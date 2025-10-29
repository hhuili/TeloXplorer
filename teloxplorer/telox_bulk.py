import os,sys
import re
import gzip
import logging
import shutil
import argparse
import itertools
import subprocess
import time
import resource
from .telox_align import extract_telomeric_reads
from . import utilities as utils
from .BloomParser import BloomParser

PRESETS = {
    'human': {
        'min_tel_freq': 0.6,
        'max_mismatch': 2,
        'max_initial_offset': 200,
        'min_tel_len': 100,
        'min_subtel_len': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': "no",
    },
    'mouse': {
        'min_tel_freq': 0.6,
        'max_mismatch': 2,
        'max_initial_offset': 200,
        'min_tel_len': 100,
        'min_subtel_len': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': "no",
    },
    'yeast': {
        'min_tel_freq': 0.7,
        'max_mismatch': 0,
        'max_initial_offset': 100,
        'min_tel_len': 30,
        'min_subtel_len': 200,
        'bloom_options': "--pwidth 10 --ydis 0.2",
        'primary_merge': "yes",
    },
    'arabidopsis': {
        'min_te_freq': 0.5,
        'max_mismatch': 2,
        'max_initial_offset': 200,
        'min_tel_len': 100,
        'min_subtel_len': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': "no",
    }
}

DEFAULTS = {
    'min_read_len': 1000,
    'min_tel_freq': 0.6,
    'max_mismatch': 2,
    'max_initial_offset': 200,
    'min_tel_len': 100,
    'min_subtel_len': 200,
    'bloom_options': "--bloom_pwidth 200 --bloom_ydis 0.4",
    'primary_merge': "no",
}


def parse_fastq_file(fastq_file,min_read_length):

    opener = gzip.open if fastq_file.endswith('.gz') else open
    mode = 'rt' if fastq_file.endswith('.gz') else 'r'

    aln_data = {}
    with opener(fastq_file, mode) as f_fastq:
        while True:
            lines = [f_fastq.readline() for _ in range(4)]
            if not lines[0]:
                break

            header, seq, plus_line, qual = lines
            read_id = header[1:].split(' ', 1)[0].strip()
            query_length = len(seq)

            if query_length < min_read_length:
                continue

            aln_data[read_id] = {'ref_name': "bulk",
                                    'ref_start': 0,
                                    'ref_end': query_length,
                                    'query_seq': seq,
                                    'strand': "+",
                                    'query_len': query_length,
                                    'query_align_len': query_length,
                                    'left_softclip': 0,
                                    'right_softclip': 0,
                                    'chrom_group': "bulkTel"}

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

def prepare_tel_segments(aln_data,tel_repeat,output_file,primary_merge):

    with open(output_file,'w') as f_out:
        for read_name in aln_data:
            chrom = aln_data[read_name]['ref_name']
            query_seq = aln_data[read_name]['query_seq']
            tel_segments, tel_repeat_freq = utils.calculate_exact_telomere_freq(query_seq,tel_repeat)
            tel_segments_update = add_nontel_segments(tel_segments,query_seq)
            if primary_merge == "yes" or sum(tel_repeat_freq) > 10000:
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

def assign_tel_group(tel_length_l,tel_length_r,has_tel_l,has_tel_r,chrom_group):

    if has_tel_l and has_tel_r:
        tel_group = f'miniTel'
        winning_tel_arm = "L,R"
        winning_tel_length = f'{tel_length_l},{tel_length_r}'
    elif has_tel_l:
        tel_group = chrom_group
        winning_tel_arm = "L"
        winning_tel_length = tel_length_l
    elif has_tel_r:
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

    bulktel_dict,minitel_dict,filtertel_dict = {},{},{}

    for read_name,aln_record in aln_data.items():
        tel_record_l = bloom_data["L"][read_name]
        tel_record_r = bloom_data["R"][read_name]
        #chrom = aln_record['ref_name']
        #mapping_start = aln_record['ref_start']
        query_align_length = aln_record['query_align_len']
        left_softclip = aln_record['left_softclip']
        right_softclip = aln_record['right_softclip']
        #mapping_end = aln_record['ref_end']
        #mapping_arm = aln_record['mapping_arm']
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

        tel_group,winning_tel_arm,winning_tel_length = assign_tel_group(tel_length_l,tel_length_r,has_tel_l,has_tel_r,chrom_group)

        final_tel_record = [read_name,
                            tel_group,
                            winning_tel_arm,
                            str(winning_tel_length)]

        if tel_group == "filterTel":
            ## filtered telomeres
            filtertel_dict[read_name] = final_tel_record
        elif tel_group == "bulkTel":
            ## bulk telomeres
            bulktel_dict[read_name] = final_tel_record
        else:
            ## mini telomeres (both end telomeres)
            minitel_dict[read_name] = final_tel_record

    winning_telomeres = {'bulktel': bulktel_dict,
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
        f_out.write(f'read_name\ttel_group\ttel_arm\ttel_length\n')
        for tel_record in tel_record_data.values():
            f_out.write('\t'.join(tel_record)+'\n')

def write_tel_seq(winning_telomeres,bloom_data,output):

    with open(output,'w') as f_out:
        for read_name, tel_record in winning_telomeres.items():
            tel_arm = tel_record[2]
            tel_length = tel_record[3]
            tel_seq = bloom_data[tel_arm][read_name].tel_seq
            f_out.write(f'>{read_name}|{tel_length}\n{tel_seq}\n')

def run_telox_bulk(fastq_input, outdir, out_prefix, args):
    """
    Alignment-free telomere length estimation.
    """

    tel_repeat_dict = utils.parse_tel_repeats(args.tel_repeats)

    bloom_data = {}
    for tel_arm in ["L", "R"]:
        tel_repeat = tel_repeat_dict[tel_arm]
        bloom_data[tel_arm] = {}

        # Paths for intermediate files
        out_tel_segments = os.path.join(outdir, f'{out_prefix}.tel_segments_{tel_arm}.txt')
        bloom_output = os.path.join(outdir, f'{out_prefix}.bloom_{tel_arm}.txt')

        read_data = parse_fastq_file(fastq_input,args.min_read_len)
        prepare_tel_segments(read_data, tel_repeat, out_tel_segments, args.primary_merge)
        run_bloom(out_tel_segments, args.bloom_options, outdir, bloom_output, tel_arm)
        seq_data = load_sequences_from_aln_reader(read_data)
        bloom_parser = BloomParser(
            tel_arm=tel_arm,
            bloom_file=bloom_output,
            seq_data=seq_data,
            tel_repeat_str=tel_repeat,
            min_labeling_threshold=args.min_tel_freq,
            max_initial_offset=args.max_initial_offset,
            fuzzy_freq_threshold=args.min_tel_freq,
            max_mismatch=args.max_mismatch
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
    winning_telomeres = get_winning_telomeres(bloom_data,read_data,args.min_tel_len,args.min_subtel_len)

    ## write final telomere reads
    tel_groups = ['bulktel','minitel','filtertel']
    for tel_group in tel_groups:
        output_telomere_length = os.path.join(outdir,f"{out_prefix}.{tel_group}_length.tsv")
        output_telomere_seq = os.path.join(outdir,f"{out_prefix}.{tel_group}_seq.fa")

        logging.info(f"> Writing {tel_group} length to: {output_telomere_length}")
        write_tel_length(winning_telomeres[tel_group],output_telomere_length)

        if tel_group in ['bulktel']:
            write_tel_seq(winning_telomeres[tel_group],bloom_data,output_telomere_seq)

def setup_common_args(parser):

    group = parser.add_argument_group("General Options")
    group.add_argument("-fq", "--long-read-fastq",
                        required=True,
                        help="Input long-read FASTQ file (gzipped or plain)")
    group.add_argument("-R", "--tel-repeats",
                       required=True,
                       help="Telomere repeat definition file (one per line: arm<TAB>repeat, regex supported)")
    group.add_argument("--out-prefix",
                       default="asm",
                       help="Output prefix [default: %(default)s].")
    group.add_argument("--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")
    group.add_argument("--start-step",
                        type=int, 
                        default=1, 
                        choices=[1, 2],
                        help="Step1: extracting and mapping telomere-like reads; "
                             "Step2: telomere length estimation [default: %(default)s]")
    group.add_argument("--debug",
                       action="store_true",
                       help="Enable debug output")

    return parser

def setup_bulk_args(parser):

    group_name = 'Telox-bulk Options'
    group = parser.add_argument_group(group_name, 'Bulk telomere length eastimation')

    group.add_argument("--preset",
                       choices=['human', 'mouse', 'yeast', 'arabidopsis'],
                       help="Species preset configuration. Manual parameters override preset values.")
    group.add_argument("-rl", "--min-read-len",
                       type=int,
                       default=None,
                       help="Minimum read length (bp) [default: %(default)s]")
    group.add_argument("-f", "--min-tel-freq",
                       type=float,
                       default=None,
                       help=f"Frequency threshold for telomere segment definition [default: {DEFAULTS['min_tel_freq']}]")
    group.add_argument("--primary-merge",
                       choices=["yes","no"],
                       default=None,
                       help=f"Merge adjacent telomere segments pre-BLOOM [default: {DEFAULTS['primary_merge']}]")
    group.add_argument("--max-mismatch",
                       type=int,
                       default=None,
                       help=f"Mismatch tolerance for telomere re-labeling [default: {DEFAULTS['max_mismatch']}]")
    group.add_argument("--max-initial-offset",
                       type=int,
                       default=None,
                       help=f"Max non-telomere length at read start [default: {DEFAULTS['max_initial_offset']}]")
    group.add_argument("-tl", "--min-tel-len",
                       type=int,
                       default=None,
                       help=f"Minimum telomere length [default: {DEFAULTS['min_tel_len']}]")
    group.add_argument("-sl", "--min-subtel-len",
                       type=int,
                       default=None,
                       help=f"Minimum sub-telomere length [default: {DEFAULTS['min_subtel_len']}]")
    group.add_argument("--bloom-options",
                       default=None,
                       help=f"BLOOM merging parameters [default: '{DEFAULTS['bloom_options']}']")

    return parser

def apply_presets_and_defaults(args):

    # 1. apply preset
    if args.preset:
        preset_values = PRESETS.get(args.preset)
        if not preset_values:
            available_presets = list(PRESETS.keys())
            raise ValueError(f"Invalid preset '{args.preset}'. Available presets are: {available_presets}")

        for key, value in preset_values.items():
            if getattr(args, key, None) is None:
                setattr(args, key, value)

    # 2. apply defaults
    for key, value in DEFAULTS.items():
        if getattr(args, key, None) is None:
            setattr(args, key, value)

    return args

def parse_arguments():

    parser = argparse.ArgumentParser(
        description="Telox-bulk: bulk telomere length eastimation"
    )

    parser = setup_common_args(parser)
    parser = setup_bulk_args(parser)

    ## apply preset and default options
    args = parser.parse_args()
    args = apply_presets_and_defaults(args)

    return args

def telox_bulk_cli():

    pipeline_start_time = time.time()

    module = "telox-bulk"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    log_file_path = os.path.join(args.outdir, f"{args.out_prefix}.log")

    utils.setup_logging(log_file_path)

    #### Log startup information
    utils.log_parameter_summary(args, module)

    #### Define output paths and file names
    paths = {
        'tel_reads_dir': "01_telomere_reads",
        'tel_length_dir': "02_telomere_length",
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    #### Intermediate files
    tel_reads_fastq = os.path.join(args.tel_reads_dir, f'{args.out_prefix}.tel_like.fastq')

    if args.start_step <= 1:
        logging.info(f"\n> STEP1: extracting telomere-like reads")
        extract_telomeric_reads(fastq_input=args.long_read_fastq,
                                tel_repeats=args.tel_repeats,
                                fastq_output=tel_reads_fastq)
    else:
        logging.info(f"\n> Skipping STEP1")

    if args.start_step <= 2:
        logging.info(f"\n> STEP2: alignment-free telomere length estimation")
        run_telox_bulk(fastq_input=tel_reads_fastq,
                       outdir=args.tel_length_dir,
                       out_prefix=args.out_prefix,
                       args=args)
    else:
        logging.info(f"\n> Skipping STEP2")

    logging.info(f"\n**** {module} completed! ****")
    total_elapsed = time.time() - pipeline_start_time

    usage = resource.getrusage(resource.RUSAGE_SELF)
    cpu_time_used = usage.ru_utime + usage.ru_stime
    logging.info(f"CPU time(s) : {cpu_time_used:.1f}")
    logging.info(f"Elapsed time(s) : {total_elapsed:.1f}")

if __name__ == "__main__":

    telox_bulk_cli()

