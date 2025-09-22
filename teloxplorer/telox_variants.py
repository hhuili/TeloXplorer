import os
import re
import argparse
import logging
from collections import defaultdict
from . import utilities as utils

def read_fasta(file_path):

    sequences = {}
    with open(file_path, 'r') as f:
        header = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    sequences[header] = "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line.upper())
        if header:
            sequences[header] = "".join(sequence)
    return sequences

def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)
    if len(s2) == 0:
        return len(s1)
    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

def get_cyclic_shifts(s):
    return [s[i:] + s[:i] for i in range(len(s))]

def get_cyclic_shifts_no_self(s):
    return [s[i:] + s[:i] for i in range(len(s)) if i != 0]

def get_canonical_form(s):
    return min(get_cyclic_shifts(s))

def find_best_representative_repeat(shifts, reference_repeat):

    min_distance = float('inf')
    best_representative = ""

    ## for inversion variants
    rc_shifts = [utils.sequence_revcomp(shift) for shift in shifts]
    all_shifts = shifts + rc_shifts

    for shift in all_shifts:
        distance = levenshtein_distance(shift, reference_repeat)

        if distance < min_distance:
            min_distance = distance
            best_representative = shift
        elif distance == min_distance:
            if best_representative == "":
                best_representative = shift
            elif shift.startswith(reference_repeat[0]):
                #best_representative = choose_longer_prefix_repeat(shift,best_representative)
                best_representative = shift
            elif shift.startswith(utils.sequence_revcomp(reference_repeat)[0]):
                #best_representative = choose_longer_prefix_repeat(shift,best_representative)
                best_representative = shift

    if best_representative in rc_shifts:
        best_representative = shifts[rc_shifts.index(best_representative)]

    return best_representative

def deduplicate_repeats(repeat_counts, reference_repeat):

    processed_canonicals = set()
    final_repeats = []

    for repeat in repeat_counts:
        canonical = get_canonical_form(repeat)
        if canonical in processed_canonicals:
            continue
        
        shifts = get_cyclic_shifts(repeat)
        best_shift = find_best_representative_repeat(shifts, reference_repeat)

        best_shifts = [best_shift]
        for _ in range(2):
            if best_shift and best_shift[0] == best_shift[-1]:
                shifts = [shift for shift in get_cyclic_shifts(best_shift) if shift not in best_shifts]
                best_shift = find_best_representative_repeat(shifts, reference_repeat)
                best_shifts.append(best_shift)

        if best_shift:
            final_repeats.append(best_shift)
        processed_canonicals.add(canonical)

    return final_repeats

def find_tvr_unit(seq, k_values, reference_repeat, min_repeat = 2):

    results = defaultdict(lambda: defaultdict(int))

    for k in k_values:
        # Iterate through 'phases' to find non-overlapping tandem repeats 
        for phase in range(k):
            last_kmer = None
            current_run_length = 0

            for start in range(phase, len(seq) - k + 1, k):
                kmer = seq[start:start+k]
                
                # Skip k-mers with non-standard bases, breaking the current run 
                if any(base not in 'ACGT' for base in kmer):
                    last_kmer = None
                    current_run_length = 0
                    continue

                if kmer == last_kmer:
                    # If the repeat is the same as the last one, extend the run
                    current_run_length += 1
                else:
                    # If the repeat changes, the run is broken; start a new one
                    last_kmer = kmer
                    current_run_length = 1
                
                # Update the overall maximum repeat length found for this canonical repeat
                if current_run_length > results[k][kmer]:
                    results[k][kmer] = current_run_length
    
    # Filter out repeats that never met the minimum repeat threshold
    filtered_results = defaultdict(dict)
    for k, repeat_counts in results.items():
        deduplicate_repeat_counts = deduplicate_repeats(repeat_counts,reference_repeat)

        for repeat in deduplicate_repeat_counts:
            if len(set(repeat)) == 1:
                continue
            count = repeat_counts[repeat]
            if repeat*min_repeat in seq:
                filtered_results[k][repeat] = count

    return filtered_results

def calculate_tvr_freq(tel_repeats,sequence,k_values,find_singleton):

    def count_singleton(sequence):
        singletons = [singleton for singleton in sequence.split('X') if singleton]

        count_dict = defaultdict(int)
        for singleton in singletons:
            count_dict[singleton] += 1

        return dict(count_dict)

    tvr = {}
    remaining_sequence = sequence
    tel_length = len(sequence)
    for k, repeats in reversed(tel_repeats.items()):
        if not repeats:
            #print(f"\nk={k}: No significant repeats found.")
            continue

        for kmer in repeats:
            kmer_matches = re.findall(kmer, remaining_sequence)
            if kmer_matches:
                counts = len(kmer_matches)
                freq = counts*k / tel_length
                tvr[kmer] = (counts,freq)

            remaining_sequence = remaining_sequence.replace(kmer, 'X'*len(kmer))

    sorted_tvr = dict(sorted(tvr.items(), key=lambda item: item[1][0], reverse=True))

    if find_singleton:
        tvr2 = {}
        singletons = count_singleton(remaining_sequence)
        for singleton, counts in singletons.items():
            singleton_len = len(singleton)
            if singleton_len < min(k_values) or singleton_len > max(k_values):
                continue
            freq = counts*len(singleton) / tel_length
            tvr2[singleton] = (counts,freq)
        sorted_tvr2 = dict(sorted(tvr2.items(), key=lambda item: item[1][0], reverse=True))
    else:
        sorted_tvr2 = {}

    return sorted_tvr, sorted_tvr2

def chrom_to_sort_key(chrom_name):
    """
    - Standard numbers: chr1, chr2, ..., chr22, chrX, chrY, chrM
    - Roman numerals: chrI, chrII, chrIII, ..., chrXVI
    - With suffixes: chr1_MATERNAL, chrII_PATERNAL, etc.
    """
    base_chrom = re.split(r'[_\W]', chrom_name)[0]

    roman_numeral = re.search(r'chr([IVXLCDM]+)', base_chrom, re.IGNORECASE)
    if roman_numeral:
        roman_str = roman_numeral.group(1).upper()
        roman_to_int = {
            'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5,
            'VI': 6, 'VII': 7, 'VIII': 8, 'IX': 9, 'X': 10,
            'XI': 11, 'XII': 12, 'XIII': 13, 'XIV': 14, 'XV': 15, 'XVI': 16
        }
        chrom_num = roman_to_int.get(roman_str, 99)
        return (1, chrom_num, chrom_name)

    standard_num = re.search(r'chr(\d+|X|Y|M)', base_chrom, re.IGNORECASE)
    if standard_num:
        num_str = standard_num.group(1)

        if num_str.upper() == 'X':
            chrom_num = 23
        elif num_str.upper() == 'Y':
            chrom_num = 24
        elif num_str.upper() == 'M':
            chrom_num = 25
        else:
            chrom_num = int(num_str)
        return (0, chrom_num, chrom_name)
    
    # other chrom format
    return (2, 0, chrom_name)

def run_telox_variants(telomere_fasta, tel_repeats, out_prefix, outdir, args):

    ####
    k_values=[int(x) for x in args.kmers.split(',')]
    min_variant_repeats=args.min_variant_repeats
    find_singleton=args.find_singleton

    tel_repeat_dict = utils.parse_tel_repeats(tel_repeats)

    records = read_fasta(telomere_fasta)
    output_variants = os.path.join(outdir,f"{out_prefix}.TVR.tsv")
    output_summary = os.path.join(outdir,f"{out_prefix}.TVR_summary.tsv")
    output_summary2 = os.path.join(outdir,f"{out_prefix}.Singleton_summary.tsv")
    if find_singleton:
        f_sum2 = open(output_summary2,'w')
        f_sum2.write(f"Chrom\tChrom_arm\tRepeat_unit\tCount\n")

    with open(output_variants, 'w') as f_out, open(output_summary, 'w') as f_sum:

        f_out.write(f"Read_name\tRepeat_unit\tRepeat_Group\tCount\tRatio\n")
        f_sum.write(f"Chrom\tChrom_arm\tRepeat_unit\tCount\n")

        tvr_summary = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        tvr_summary2 = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

        for read_name, sequence in records.items():
            tel_arm = read_name.split('|')[-2]
            chrom = read_name.split('|')[-3]
            canonical_repeat = tel_repeat_dict[tel_arm]

            tel_repeat_units = find_tvr_unit(sequence, k_values, canonical_repeat, min_variant_repeats)
            tvr_per_read, tvr_per_read2 = calculate_tvr_freq(tel_repeat_units, sequence, k_values, find_singleton)

            ## write TVR per read
            for repeat_unit, (count, freq) in tvr_per_read.items():
                if count > 1:
                    tag = "non_singleton"
                    f_out.write(f'{read_name}\t{repeat_unit}\t{tag}\t{count}\t{freq}\n')
                    tvr_summary[chrom][tel_arm][repeat_unit] += count
                else:
                    tag = "singleton"
                    f_out.write(f'{read_name}\t{repeat_unit}\t{tag}\t{count}\t{freq}\n')
                    tvr_summary2[chrom][tel_arm][repeat_unit] += count

            # singleton
            if find_singleton and tvr_per_read2:
                for repeat_unit, (count, freq) in tvr_per_read2.items():
                    tag = "singleton"
                    f_out.write(f'{read_name}\t{repeat_unit}\t{tag}\t{count}\t{freq}\n')
                    tvr_summary2[chrom][tel_arm][repeat_unit] += count

        for chrom in sorted(tvr_summary.keys(), key=chrom_to_sort_key):
            for tel_arm in sorted(tvr_summary[chrom].keys()):
                sorted_repeats = sorted(
                    tvr_summary[chrom][tel_arm].items(),
                    key=lambda x: (-x[1], x[0])
                )

                for repeat_unit, count in sorted_repeats:
                    f_sum.write(f"{chrom}\t{tel_arm}\t{repeat_unit}\t{count}\n")

        if find_singleton:
            for chrom in sorted(tvr_summary2.keys(), key=chrom_to_sort_key):
                for tel_arm in sorted(tvr_summary2[chrom].keys()):
                    sorted_repeats = sorted(
                        tvr_summary2[chrom][tel_arm].items(),
                        key=lambda x: (-x[1], x[0])
                    )

                    for repeat_unit, count in sorted_repeats:
                        f_sum2.write(f"{chrom}\t{tel_arm}\t{repeat_unit}\t{count}\n")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='telox_variants: Teloxplorer module for chromosome-specific telomere variant repeat analysis.'
    )
    parser.add_argument(
        "-i", "--telomere-fasta",
        required=True,
        help='Input telomere sequence file',
    )
    parser.add_argument(
        "-r", "--tel-repeats",
        required=True,
        help="Path to telomere repeat definition file (one repeat per line, format: arm<TAB>sequence, supports regex).",
    )
    parser.add_argument(
        "-k", "--kmers",
        type=str,
        default="5,6,7",
        help="Kmer sizes to find telomeric repeat units (default: 5,6,7).",
    )
    parser.add_argument(
        "--min-variant-repeats",
        type=int,
        default=2,
        help="Minimum number of consecutive repeats for telomere variant repeats calling (default: 2).",
    )
    parser.add_argument(
        "--find-singleton",
        action='store_true',
        required=False,
        help="Find singleton telomeric repeat units.",
    )
    parser.add_argument(
        "-o", "--out-prefix",
        default="telox",
        help="Output prefix (default: telox).",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory (default: .).",
    )

    return parser.parse_args()

def variants_cli():
    """
    The main entry point for the 'telox-variants' command-line tool.
    """

    utils.setup_logging()
    args = parse_arguments()

    #### Define output paths and file names
    paths = {
        'tel_variants_dir': "telomere_variants"
    }

    for name, subdir in paths.items():
        path = os.path.join(args.outdir, subdir)
        args = utils.add_path(args, name, path)

    logging.info(f"Analyzing telomere variant repeats.")
    run_telox_variants(telomere_fasta=args.telomere_fasta,
                    tel_repeats=args.tel_repeats,
                    out_prefix=args.out_prefix,
                    outdir=args.tel_variants_dir,
                    args=args)


if __name__ == "__main__":
    variants_cli()
