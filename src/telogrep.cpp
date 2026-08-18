/*
   Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>
   Licensed under the GNU General Public License v3.0
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <zlib.h>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <memory>
#include <re2/re2.h>
#include <cmath>
#include <cstdio>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

struct FastqRecord {
    std::string name;
    std::string comment;
    std::string sequence;
    std::string quality;
};

std::string repeat_string(const std::string& s, int n) {
    if (n <= 0) return "";
    std::string result;
    result.reserve(s.length() * n);
    for (int i = 0; i < n; ++i) {
        result += s;
    }
    return result;
}

void print_usage() {
    std::cerr << "\n";
    std::cerr << "Usage:   telogrep [options]\n\n";
    std::cerr << "Input/Output (Defaults to stdin/stdout):\n";
    std::cerr << "  -i, --input <file>       Input FASTQ file [stdin].\n";
    std::cerr << "  -o, --output <file>      Output FASTQ file [stdout].\n\n";
    std::cerr << "Filtering:\n";
    std::cerr << "  -l, --min-len <int>      Discard reads shorter than INT [0].\n";
    std::cerr << "  -q, --min-qual <float>   Discard reads with avg quality below FLOAT [0.0].\n\n";
    std::cerr << "Motif Search:\n";
    std::cerr << "  -L, --left <str>         Left arm motif UNIT (e.g., 'CCCTAA' or regex 'C{1,3}A').\n";
    std::cerr << "  -R, --right <str>        Right arm motif UNIT (e.g., 'TTAGGG' or regex 'TG{1,3}').\n";
    std::cerr << "  -n, --min-copy <int>     Minimum number of motif copies required [1]. Applies to both exact and regex modes.\n";
    std::cerr << "  --regex                  Treat L-/R-motif as regex. Program builds '(motif){n,}' automatically.\n";
    std::cerr << "General:\n";
    std::cerr << "  -t, --threads <int>      Number of worker threads [4].\n";
    std::cerr << "  -h, --help               Show this help message.\n\n";
}

int main(int argc, char* argv[]) {
    if (argc == 1) { print_usage(); return 1; }

    std::string input_path, output_path, pattern_L, pattern_R;
    int min_copy = 1;
    int num_threads = 4;
    bool regex_mode = false;
    size_t min_len = 0;
    double min_qual = 0.0;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];

            auto next_arg = [&](const char* opt_name) {
                if (++i >= argc) throw std::runtime_error(std::string(opt_name) + " requires a value.");
                return std::string(argv[i]);
            };

            if (arg == "--input" || arg == "-i") input_path = next_arg("--input/-i");
            else if (arg == "--output" || arg == "-o") output_path = next_arg("--output/-o");
            else if (arg == "--left" || arg == "-L") pattern_L = next_arg("--left/-L");
            else if (arg == "--right" || arg == "-R") pattern_R = next_arg("--right/-R");
            else if (arg == "--min-copy" || arg == "-n") min_copy = std::stoi(next_arg("--min-copy/-n"));
            else if (arg == "--threads" || arg == "-t") num_threads = std::stoi(next_arg("--threads/-t"));
            else if (arg == "--min-len" || arg == "-l") min_len = std::stoul(next_arg("--min-len/-l"));
            else if (arg == "--min-qual" || arg == "-q") min_qual = std::stod(next_arg("--min-qual/-q"));
            else if (arg == "--regex") regex_mode = true;
            else if (arg == "--help" || arg == "-h") { print_usage(); return 0; }
            else { std::cerr << "[Error] Unknown argument: " << arg << std::endl; print_usage(); return 1; }
        }
    } catch (const std::exception& e) {
        std::cerr << "Argument Parsing Error: " << e.what() << std::endl;
        print_usage();
        return 1;
    }

    bool has_LR = !pattern_L.empty() && !pattern_R.empty();
    bool filter_active = min_len > 0 || min_qual > 0;

    bool motif_search_active = has_LR;

    if (motif_search_active) {
        if (pattern_L.empty() || pattern_R.empty()) {
            std::cerr << "Error: Both --left (-L) and --right (-R) motifs are required for search.\n";
            return 1;
        }
        if (min_copy < 1) {
            std::cerr << "Error: --min-copy (-n) must be >= 1.\n";
            return 1;
        }
    }

    if (!motif_search_active && !filter_active) {
        std::cerr << "Error: No action specified. Please set a filter (--min-len/--min-qual) or a motif search.\n";
        print_usage();
        return 1;
    }

    if(filter_active) {
        std::cerr << "Filter mode: ON";
        if (min_len > 0) std::cerr << " | min-len: " << min_len;
        if (min_qual > 0) std::cerr << " | min-qual: " << std::fixed << std::setprecision(1) << min_qual;
        std::cerr << std::endl;
    }
    if(motif_search_active) {
         std::cerr << "Motif search mode: ON | Regex: " << (regex_mode ? "ON" : "OFF") << std::endl;
    }

    auto start_time = std::chrono::high_resolution_clock::now();

    try {
        // Pre-compute Phred lookup table
        std::vector<double> phred_lookup(128);
        for (int i = 0; i < 128; ++i) {
            phred_lookup[i] = std::pow(10.0, -static_cast<double>(i) / 10.0);
        }

        std::string search_str_L, search_str_R;
        std::unique_ptr<re2::RE2> regex_L, regex_R;

        if (motif_search_active) {
            if (!regex_mode) {
                // Exact match mode: construct the long string
                search_str_L = repeat_string(pattern_L, min_copy);
                search_str_R = repeat_string(pattern_R, min_copy);
            } else {
                // Regex mode: construct the regex (Unit){n,}
                search_str_L = "(" + pattern_L + "){" + std::to_string(min_copy) + ",}";
                search_str_R = "(" + pattern_R + "){" + std::to_string(min_copy) + ",}";

                regex_L = std::make_unique<re2::RE2>(search_str_L);
                regex_R = std::make_unique<re2::RE2>(search_str_R);

                if (!regex_L->ok()) {
                    throw std::runtime_error("Invalid Left Regex: " + regex_L->error() + "\nGenerated: " + search_str_L);
                }
                if (!regex_R->ok()) {
                    throw std::runtime_error("Invalid Right Regex: " + regex_R->error() + "\nGenerated: " + search_str_R);
                }
            }
        }

        gzFile input_file;
        if (input_path.empty() || input_path == "-") {
            input_file = gzdopen(fileno(stdin), "r");
        } else {
            input_file = gzopen(input_path.c_str(), "rb");
        }
        if (!input_file) throw std::runtime_error("Could not open input file (or stdin)");

        std::ofstream file_out;
        std::ostream* out_stream_ptr = nullptr;
        if (output_path.empty() || output_path == "-") {
            out_stream_ptr = &std::cout;
        } else {
            file_out.open(output_path);
            if (!file_out.is_open()) {
                gzclose(input_file);
                throw std::runtime_error("Could not open output file");
            }
            out_stream_ptr = &file_out;
        }

        kseq_t *seq = kseq_init(input_file);

        omp_set_num_threads(num_threads);
        
        long long total_reads_found = 0;
        const size_t BATCH_SIZE_READS = 5000;
        bool eof = false;

        std::vector<FastqRecord> records_batch;
        records_batch.reserve(BATCH_SIZE_READS);

        std::vector<std::string> output_buffers(num_threads);
        for(auto& buf : output_buffers) {
            buf.reserve((BATCH_SIZE_READS * 1000) / num_threads);
        }

        while (!eof) {
            records_batch.clear(); 

            for (size_t i = 0; i < BATCH_SIZE_READS; ++i) {
                int status = kseq_read(seq);
                if (status < 0) {
                    eof = true;
                    break;
                }
                if (!seq->seq.l || !seq->qual.l) continue;

                records_batch.push_back({
                    seq->name.l > 0 ? std::string(seq->name.s, seq->name.l) : "",
                    seq->comment.l > 0 ? std::string(seq->comment.s, seq->comment.l) : "",
                    seq->seq.l > 0 ? std::string(seq->seq.s, seq->seq.l) : "",
                    seq->qual.l > 0 ? std::string(seq->qual.s, seq->qual.l) : ""
                });
            }

            if (records_batch.empty()) break;

            for(auto& buf : output_buffers) buf.clear();
            
            long long reads_found_batch = 0;

            #pragma omp parallel reduction(+:reads_found_batch)
            {
                int thread_id = omp_get_thread_num();
                std::string& local_buffer = output_buffers[thread_id]; 
                
                #pragma omp for schedule(dynamic, 100)
                for (size_t i = 0; i < records_batch.size(); ++i) {
                    const auto& record = records_batch[i];
                    
                    bool passed_filters = true;
                    bool found_motif = false;

                    // Length & Quality Filters
                    if (filter_active) {
                        if (min_len > 0 && record.sequence.length() < min_len) {
                            passed_filters = false;
                        }

                        if (passed_filters && min_qual > 0) {
                            double error_prob_sum = 0.0;
                            for (char c : record.quality) {
                                int phred_score = c - 33;
                                if (phred_score < 0) phred_score = 0;
                                if (phred_score > 127) phred_score = 127; 
                                error_prob_sum += phred_lookup[phred_score];
                            }
                            double avg_error_prob = error_prob_sum / record.sequence.length();
                            if (avg_error_prob > 0.0) {
                                double avg_qual_score = -10.0 * std::log10(avg_error_prob);
                                if (avg_qual_score < min_qual) passed_filters = false;
                            }
                        }
                    }

                    if (filter_active && !passed_filters) continue;

                    // Motif Search
                    if (motif_search_active) {
                        if (regex_mode) {
                            if (re2::RE2::PartialMatch(record.sequence, *regex_L) || re2::RE2::PartialMatch(record.sequence, *regex_R)) {
                                found_motif = true;
                            }
                        } else {
                            if (record.sequence.find(search_str_L) != std::string::npos || record.sequence.find(search_str_R) != std::string::npos) {
                                found_motif = true;
                            }
                        }
                    }

                    bool write_this_read = (!motif_search_active && passed_filters) || (motif_search_active && found_motif);
                    
                    if (write_this_read) {
                        local_buffer += '@';
                        local_buffer += record.name;
                        if (!record.comment.empty()) {
                            local_buffer += ' ';
                            local_buffer += record.comment;
                        }
                        local_buffer += '\n';
                        local_buffer += record.sequence;
                        local_buffer += "\n+\n";
                        local_buffer += record.quality;
                        local_buffer += '\n';
                        
                        reads_found_batch++;
                    }
                }
            }

            for (const auto& local_buffer : output_buffers) {
                if (!local_buffer.empty()) {
                    out_stream_ptr->write(local_buffer.c_str(), local_buffer.length());
                }
            }
            total_reads_found += reads_found_batch;
        }

        kseq_destroy(seq);
        if (!input_path.empty() && input_path != "-") gzclose(input_file);
        if (!output_path.empty() && output_path != "-") file_out.close();

        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end_time - start_time;
        
        //std::cerr << "Finished in " << std::fixed << std::setprecision(2) << duration.count() << " seconds." << std::endl;
        std::cerr << "Found and wrote " << total_reads_found << " reads." << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
