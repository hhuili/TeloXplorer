import os
import logging
from .utilities import run_cmd

class TelomereMapper:
    """
    A class to encapsulate the telomere reads mapping pipeline.

    This class runs minimap2 to align reads, then uses samtools to sort,
    index, and generate statistics for the resulting alignment.
    """

    def __init__(self, fastq_input, ref_genome, out_prefix, out_dir,
                 threads, mm2_preset, min_mapping_quality):
        """
        Initializes the TelomereMapper pipeline.

        Args:
            out_prefix (str): The unique identifier for the sample.
            fastq_input (str): Path to the FASTQ file.
            ref_genome (str): Path to the reference genome FASTA file.
            out_dir (str): Directory to save output files.
            threads (int, optional): Number of threads to use.
            mm2_preset (str, optional): minimap2 preset options.
            min_mapping_quality (int, optional): Minimum mapping quality for samtools.
        """
        # Store configuration parameters
        self.fastq_input = fastq_input
        self.ref_genome = ref_genome
        self.out_prefix = out_prefix
        self.out_dir = out_dir
        self.threads = str(threads)
        self.mm2_preset = mm2_preset
        self.min_mapping_quality = str(min_mapping_quality)
        self.max_secondary_alignments = '2'

        # Executable paths (can be modified if not in system PATH)
        self.minimap2_exec = 'minimap2'
        self.samtools_exec = 'samtools'

        # Define all output file paths
        base_path = os.path.join(self.out_dir, self.out_prefix)
        self.bam_raw_path = f"{base_path}.bam"
        self.bam_sorted_path = f"{base_path}.sort.bam"
        self.stats_path = f"{base_path}.samstat"

        # Ensure output directory exists
        os.makedirs(self.out_dir, exist_ok=True)

    def _run_mapping(self):
        """Step 1: Run minimap2 alignment and samtools view."""
        command = [
            self.minimap2_exec,
            '-t', self.threads,
            self.mm2_preset,
            '-N', self.max_secondary_alignments,
            '-Y',
            self.ref_genome,
            self.fastq_input,
            '|',
            self.samtools_exec, 'view',
            '-bS',
            '-q', self.min_mapping_quality,
            '-o', self.bam_raw_path
        ]
        run_cmd(command, use_shell=True)

    def _sort_bam(self):
        """Step 2: Sort the generated BAM file."""
        command = [self.samtools_exec, 'sort', '-o', self.bam_sorted_path, self.bam_raw_path]
        run_cmd(command)

    def _index_bam(self):
        """Step 3: Index the sorted BAM file."""
        command = [self.samtools_exec, 'index', self.bam_sorted_path]
        run_cmd(command)

    def _generate_stats(self):
        """Step 4: Generate alignment statistics."""
        command = [self.samtools_exec,'flagstat',self.bam_sorted_path,'>',self.stats_path]
        run_cmd(command, use_shell=True)

    def _cleanup(self):
        """Step 5: Remove the initial, unsorted BAM file."""
        logging.info(f"> Cleaning up intermediate file: {self.bam_raw_path}")
        os.remove(self.bam_raw_path)

    def run(self):
        """
        Executes the full mapping pipeline from start to finish.
        """
        self._run_mapping()
        self._sort_bam()
        self._index_bam()
        self._generate_stats()
        self._cleanup()


