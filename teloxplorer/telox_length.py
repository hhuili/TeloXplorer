import os,sys
import time
import resource
import logging
import shutil
import pysam
import argparse
import itertools
import subprocess
from . import LengthPlotter
from . import utilities as utils
from .BloomParser import BloomParser

PRESETS = {
    'human': {
        'genome_subtel_range': 500000,
        'min_tel_freq': 0.6,
        'max_mismatch': 2,
        'max_initial_offset': 200,
        'min_tel_len': 100,
        'min_subtel_len': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': "no",
        'kmers': "5,6,7"
    },
    'mouse': {
        'genome_subtel_range': 500000,
        'min_tel_freq': 0.6,
        'max_mismatch': 2,
        'max_initial_offset': 200,
        'min_tel_len': 100,
        'min_subtellen': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': "no",
        'kmers': "5,6,7"
    },
    'yeast': {
        'genome_subtel_range': 20000,
        'min_tel_freq': 0.7,
        'max_mismatch': 0,
        'max_initial_offset': 100,
        'min_tel_length': 30,
        'min_subtel_length': 200,
        'bloom_options': "--pwidth 10 --ydis 0.2",
        'primary_merge': "yes",
        'kmers': "3,4"
    },
    'arabidopsis': {
        'genome_subtel_range': 20000,
        'min_tel_freq': 0.5,
        'max_mismatch': 2,
        'max_initial_offset': 200,
        'min_tel_len': 100,
        'min_subtel_len': 200,
        'bloom_options': "--pwidth 200 --ydis 0.4",
        'primary_merge': "no",
        'kmers': "6,7,8"
    }
}

DEFAULTS = {
    'genome_subtel_range': 500000,
    'min_read_len': 1000,
    'min_tel_freq': 0.6,
    'max_mismatch': 2,
    'max_initial_offset': 200,
    'min_tel_len': 100,
    'min_subtel_len': 200,
    'bloom_options': "--bloom_pwidth 200 --bloom_ydis 0.4",
    'primary_merge': "no",
}

def add_mapping_arm(ref_name, ref_start, ref_end, chrom_sizes, genome_subtelomere_region):

    chrom_size = chrom_sizes[ref_name]

    if ref_name not in chrom_sizes:
        chrom_group = "unknown_chrom"
        mapping_arm = "unknown_chrom"
        return

    if ref_start <= genome_subtelomere_region:
        chrom_group = "chromTerminal"
        mapping_arm = "L"
    elif ref_end >= (chrom_size - genome_subtelomere_region):
        chrom_group = "chromTerminal"
        mapping_arm = "R"
    else:
        chrom_group = "chromInternal"
        mapping_arm = "none"

    return chrom_group, mapping_arm

def parse_alignment_file(bam_file,min_read_length,chrom_sizes,genome_subtelomere_region):

    with pysam.AlignmentFile(bam_file, "rb") as bamfile:
        aln_data = {}
        for aln in bamfile.fetch(until_eof=True):
            if (aln.is_unmapped or
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
                                                       genome_subtelomere_region)

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
                #logging.info(f'> Duplicated read: {read_id}')
                pass

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
    """Write telomere length to file."""

    with open(output,'w') as f_out:
        f_out.write(f'read_name\ttel_group\tchrom\ttel_arm\ttel_length\tmapping_start\tmapping_end\n')

        for tel_record in tel_record_data.values():
            f_out.write('\t'.join(tel_record)+'\n')

def write_tel_summary(tel_record_data, output, out_prefix):
    """Write telomere length summary statistics to file."""

    groups = {}
    all_lengths = []
    paternal_lengths = []
    maternal_lengths = []
    tel_group = None
    
    for record in tel_record_data.values():
        if tel_group is None:
            tel_group = record[1]
        
        chrom, arm, length = record[2], record[3], int(record[4])
        key = (chrom, arm)
        
        if key not in groups:
            groups[key] = []
        groups[key].append(length)
        all_lengths.append(length)

        if '_' in chrom:
            hap_type = chrom.split('_')[-1]
            if hap_type.upper() in ['PATERNAL', 'PAT', 'P']:
                paternal_lengths.append(length)
            elif hap_type.upper() in ['MATERNAL', 'MAT', 'M']:
                maternal_lengths.append(length)

    with open(output, 'w') as f_out:
        f_out.write("level\tsample\tchrom\ttel_arm\ttel_count\ttel_mean\ttel_median\n")

        if all_lengths:
            total_count = len(all_lengths)
            total_mean = sum(all_lengths) / total_count
            sorted_all = sorted(all_lengths)
            total_median = (sorted_all[total_count//2] if total_count % 2 == 1 
                          else (sorted_all[total_count//2 - 1] + sorted_all[total_count//2]) / 2)
            
            f_out.write(f"total\t{out_prefix}\tall\tall\t{total_count}\t{total_mean:.0f}\t{total_median:.0f}\n")

        if paternal_lengths:
            pat_count = len(paternal_lengths)
            pat_mean = sum(paternal_lengths) / pat_count
            sorted_pat = sorted(paternal_lengths)
            pat_median = (sorted_pat[pat_count//2] if pat_count % 2 == 1 
                         else (sorted_pat[pat_count//2 - 1] + sorted_pat[pat_count//2]) / 2)
            
            f_out.write(f"total\t{out_prefix}\tpaternal\tall\t{pat_count}\t{pat_mean:.0f}\t{pat_median:.0f}\n")
        
        if maternal_lengths:
            mat_count = len(maternal_lengths)
            mat_mean = sum(maternal_lengths) / mat_count
            sorted_mat = sorted(maternal_lengths)
            mat_median = (sorted_mat[mat_count//2] if mat_count % 2 == 1 
                         else (sorted_mat[mat_count//2 - 1] + sorted_mat[mat_count//2]) / 2)
            
            f_out.write(f"total\t{out_prefix}\tmaternal\tall\t{mat_count}\t{mat_mean:.0f}\t{mat_median:.0f}\n")
        
        # sort chromosomes
        sorted_chroms = sorted(groups.keys(), key=lambda x: utils.natural_sort_key(x[0]))
        
        for (chrom, arm) in sorted_chroms:
            lengths = groups[(chrom, arm)]
            count = len(lengths)
            mean_len = sum(lengths) / count if count > 0 else 0
            sorted_lengths = sorted(lengths)
            
            if count % 2 == 1:
                median_len = sorted_lengths[count // 2]
            else:
                median_len = (sorted_lengths[count//2 - 1] + sorted_lengths[count//2]) / 2
            
            f_out.write(f"chromosome\t{out_prefix}\t{chrom}\t{arm}\t{count}\t{mean_len:.0f}\t{median_len:.0f}\n")


def write_tel_seq(winning_telomeres,bloom_data,output):
    """Write telomere sequence to file."""

    with open(output,'w') as f_out:
        for read_name, tel_record in winning_telomeres.items():
            chrom = tel_record[2]
            tel_arm = tel_record[3]
            tel_length = tel_record[4]
            tel_seq = bloom_data[tel_arm][read_name].tel_seq
            f_out.write(f'>{read_name}|{chrom}|{tel_arm}|{tel_length}\n{tel_seq}\n')

def run(bam_file, outdir, out_prefix, args):
    """Alignment-based telomere length estimation."""

    tel_repeat_dict = utils.parse_tel_repeats(args.tel_repeats)
    chrom_sizes = utils.parse_genome_indx(f"{args.ref_genome}.fai")

    bloom_data = {}
    for tel_arm in ["L", "R"]:
        tel_repeat = tel_repeat_dict[tel_arm]
        bloom_data[tel_arm] = {}

        # intermediate files
        out_tel_segments = os.path.join(outdir, f'{out_prefix}.tel_segments_{tel_arm}.txt')
        bloom_output = os.path.join(outdir, f'{out_prefix}.bloom_{tel_arm}.txt')

        aln_data = parse_alignment_file(bam_file,args.min_read_len,chrom_sizes,args.genome_subtel_range)
        prepare_tel_segments(aln_data, tel_repeat, out_tel_segments, args.primary_merge)
        run_bloom(out_tel_segments, args.bloom_options, outdir, bloom_output, tel_arm)
        #segment_generator = prepare_tel_segments_generator(aln_data, tel_repeat)
        #bloom_output_str = run_bloom_in_memory(segment_generator, bloom_options)
        seq_data = load_sequences_from_aln_reader(aln_data)
        #bloom_output_stream = io.StringIO(bloom_output_str)
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
        if not args.debug:
            if os.path.exists(out_tel_segments):
                os.remove(out_tel_segments)
            os.remove(bloom_output)
        shutil.rmtree(f'{outdir}/bloom_{tel_arm}')

    ## calculate telomere length and filtering
    winning_telomeres = get_winning_telomeres(bloom_data,aln_data,args.min_tel_len,args.min_subtel_len)

    ## write final telomere reads
    tel_groups = ['chromtel','neotel','minitel','filtertel']
    for tel_group in tel_groups:
        output_telomere_length = os.path.join(outdir,f"{out_prefix}.{tel_group}_length.tsv")
        output_telomere_seq = os.path.join(outdir,f"{out_prefix}.{tel_group}_seq.fasta")

        logging.info(f"> Writing {tel_group} length to: {output_telomere_length}")
        write_tel_length(winning_telomeres[tel_group],output_telomere_length)

        if tel_group in ['chromtel','neotel']:
            write_tel_seq(winning_telomeres[tel_group],bloom_data,output_telomere_seq)

    for tel_group in ['chromtel']:
        output_telomere_summary = os.path.join(outdir,f"{out_prefix}.{tel_group}_summary.tsv")
        output_telomere_length = os.path.join(outdir,f"{out_prefix}.{tel_group}_length.tsv")

        logging.info(f"> Writing {tel_group} summary to: {output_telomere_summary}")
        write_tel_summary(winning_telomeres[tel_group],output_telomere_summary,out_prefix)

        if args.plot:
            logging.info(f"> Plotting {tel_group} telomere length")
            LengthPlotter.run_length_plot(
                length_input=output_telomere_length,
                outdir=outdir,
                out_prefix=f'{out_prefix}.{tel_group}',
                args=args
            )


def setup_common_args(parser):

    group = parser.add_argument_group("General Options")
    group.add_argument("-b", "--bam",
                       required=True,
                       help="Minimap2 aligned BAM/SAM file")
    group.add_argument("-r", "--ref-genome",
                       required=True,
                       help="Reference genome file")
    group.add_argument("-R", "--tel-repeats",
                       required=True,
                       help="Telomere repeat definition file (one per line: arm<TAB>repeat, regex supported)")
    group.add_argument("--out-prefix",
                       default="asm",
                       help="Output prefix [default: %(default)s].")
    group.add_argument("--outdir",
                       default=".",
                       help="Output directory [default: %(default)s].")
    group.add_argument("--debug",
                       action="store_true",
                       help="Enable debug output")

    return parser

def setup_length_args(parser):

    group_name = 'Telox-length Options'
    group = parser.add_argument_group(group_name, 'Chromosome-specific telomere length estimation')

    group.add_argument("--preset",
                       choices=['human', 'mouse', 'yeast', 'arabidopsis'],
                       help="Species preset configuration. Manual parameters override preset values.")
    group.add_argument("-B", "--genome-subtel-range",
                        type=int,
                        default=None,
                        help=f"Size of the subtelomeric region (bp) from the chromosome ends. "
                             f"Telomeres within this range are classified as 'terminal', "
                             f"while those outside are 'internal'. "
                             f"[default: {DEFAULTS['genome_subtel_range']}]")
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
    group.add_argument("--plot",
                        action='store_true',
                        help="Generate telomere length plot [default: False]")

    parser = LengthPlotter.setup_length_plotting_args(parser)

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
        description="Telox-length: chromosome-specific telomere length analysis"
    )

    parser = setup_common_args(parser)
    parser = setup_length_args(parser)

    ## apply preset and default options
    args = parser.parse_args()
    args = apply_presets_and_defaults(args)

    return args


def telox_length_cli():

    pipeline_start_time = time.time()

    module="telox-length"
    args = parse_arguments()

    os.makedirs(args.outdir, exist_ok=True)
    log_file_path = os.path.join(args.outdir, f"{args.out_prefix}.log")

    utils.setup_logging(log_file_path)

    #### Log startup information
    utils.log_parameter_summary(args, module)

    logging.info(f"\n> Alignment-based telomere length estimation")
    run(bam_file=args.bam,
        outdir=args.outdir,
        out_prefix=args.out_prefix,
        args=args)

    logging.info(f"\n**** {module} completed! ****")

    total_elapsed = time.time() - pipeline_start_time
    usage = resource.getrusage(resource.RUSAGE_SELF)
    cpu_time_used = usage.ru_utime + usage.ru_stime

    logging.info(f"CPU time(s) : {cpu_time_used:.1f}")
    logging.info(f"Elapsed time(s) : {total_elapsed:.1f}")

if __name__ == "__main__":

    telox_length_cli()

