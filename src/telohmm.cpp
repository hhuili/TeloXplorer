/*
   Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>
   Licensed under the GNU General Public License v3.0
*/

#include "telohmm.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

const float LOG_ZERO = -1e10f;

void ViterbiWorkspace::ensure_capacity(int N, int L) {
    L_plus_1 = L + 1;

    size_t required_size = static_cast<size_t>(N + 1) * static_cast<size_t>(L_plus_1);
    
    if (required_size > V_M.size()) {
        V_M.resize(required_size); V_I.resize(required_size); V_D.resize(required_size);
        ptr_M.resize(required_size); ptr_I.resize(required_size); ptr_D.resize(required_size);
    }

    for (int i = 0; i <= N; ++i) {
        V_M[i * L_plus_1] = LOG_ZERO;
        V_I[i * L_plus_1] = LOG_ZERO;
        V_D[i * L_plus_1] = LOG_ZERO;
    }
    for (int j = 0; j <= L; ++j) {
        V_M[j] = LOG_ZERO;
        V_I[j] = LOG_ZERO;
        V_D[j] = LOG_ZERO;
    }
    V_M[0] = 0.0f; // V_M[0][0]
}

ProfileHMM::ProfileHMM(const std::string& backbone, const std::vector<char>& valid_chars, char ref_char) {
    this->ref_char = ref_char;
    this->backbone = backbone;
    L = backbone.length();
    alphabet = valid_chars;
    alpha_size = alphabet.size();

    char_to_idx_table.assign(256, 0);
    for (size_t i = 0; i < alphabet.size(); ++i) {
        char_to_idx_table[static_cast<unsigned char>(alphabet[i])] = i;
    }
    
    allocate_matrices();
    initialize_probabilities(backbone);
}

void ProfileHMM::allocate_matrices() {
    t_MM.assign(L + 1, LOG_ZERO); t_MI.assign(L + 1, LOG_ZERO); t_MD.assign(L + 1, LOG_ZERO);
    t_IM.assign(L + 1, LOG_ZERO); t_II.assign(L + 1, LOG_ZERO);
    t_DM.assign(L + 1, LOG_ZERO); t_DD.assign(L + 1, LOG_ZERO);
    
    e_M.assign((L + 1) * alpha_size, LOG_ZERO);
    e_I.assign((L + 1) * alpha_size, LOG_ZERO);
}

void ProfileHMM::initialize_probabilities(const std::string& backbone) {
    float p_MM = 0.75f, p_MI = 0.10f, p_MD = 0.15f; 
    float p_IM = 0.70f, p_II = 0.30f;
    float p_DM = 0.50f, p_DD = 0.50f;

    for (int j = 0; j <= L; ++j) {
        t_MM[j] = std::log(p_MM); t_MI[j] = std::log(p_MI); t_MD[j] = std::log(p_MD);
        t_IM[j] = std::log(p_IM); t_II[j] = std::log(p_II);
        t_DM[j] = std::log(p_DM); t_DD[j] = std::log(p_DD);
    }

    float match_bonus = 0.50f;
    float mismatch_penalty = (alpha_size > 1) ? (1.0f - match_bonus) / (alpha_size - 1) : 1.0f;

    for (int j = 1; j <= L; ++j) {
        char expected_char = backbone[j - 1];
        int expected_idx = get_char_idx(expected_char);

        for (int k = 0; k < alpha_size; ++k) {
            int idx = j * alpha_size + k;
            if (k == expected_idx) {
                e_M[idx] = std::log(match_bonus); 
            } else {
                e_M[idx] = std::log(mismatch_penalty);
            }
            e_I[idx] = std::log(1.0f / alpha_size); 
        }
    }
}

AlignmentPath ProfileHMM::viterbi_align(const std::string& read, ViterbiWorkspace& ws) const {
    int N = read.length();
    ws.ensure_capacity(N, L);

    for (int j = 1; j <= L; ++j) {
        float score_MD = ws.get_VM(0, j-1) + t_MD[j-1];
        float score_DD = ws.get_VD(0, j-1) + t_DD[j-1];
        ws.get_VD(0, j) = std::max(score_MD, score_DD);
        ws.get_ptrD(0, j) = (score_MD > score_DD) ? 0 : 2; 
    }

    for (int i = 1; i <= N; ++i) {
        int char_idx = get_char_idx(read[i - 1]);

        for (int j = 1; j <= L; ++j) {
            float score_MM = ws.get_VM(i-1, j-1) + t_MM[j-1];
            float score_IM = ws.get_VI(i-1, j-1) + t_IM[j-1];
            float score_DM = ws.get_VD(i-1, j-1) + t_DM[j-1];
            float max_M = std::max({score_MM, score_IM, score_DM});
            
            ws.get_VM(i, j) = e_M[j * alpha_size + char_idx] + max_M;
            
            if (max_M == score_MM) ws.get_ptrM(i, j) = 0;
            else if (max_M == score_IM) ws.get_ptrM(i, j) = 1;
            else ws.get_ptrM(i, j) = 2;

            float score_MI = ws.get_VM(i-1, j) + t_MI[j];
            float score_II = ws.get_VI(i-1, j) + t_II[j];
            float max_I = std::max(score_MI, score_II);
            
            ws.get_VI(i, j) = e_I[j * alpha_size + char_idx] + max_I;
            ws.get_ptrI(i, j) = (max_I == score_MI) ? 0 : 1;

            float score_MD = ws.get_VM(i, j-1) + t_MD[j-1];
            float score_DD = ws.get_VD(i, j-1) + t_DD[j-1];
            
            ws.get_VD(i, j) = std::max(score_MD, score_DD);
            ws.get_ptrD(i, j) = (score_MD > score_DD) ? 0 : 2;
        }
    }

    float max_final_score = LOG_ZERO;
    int best_j = L;
    State current_state = M;

    for (int j = 1; j <= L; ++j) {
        if (ws.get_VM(N, j) >= max_final_score) { 
            max_final_score = ws.get_VM(N, j); best_j = j; current_state = M; 
        }
        if (ws.get_VI(N, j) >= max_final_score) { 
            max_final_score = ws.get_VI(N, j); best_j = j; current_state = I; 
        }
        if (ws.get_VD(N, j) >= max_final_score) { 
            max_final_score = ws.get_VD(N, j); best_j = j; current_state = D; 
        }
    }

    std::vector<State> reversed_path;
    reversed_path.reserve(N + L);
    int i = N, j = best_j;

    while (i > 0 || j > 0) {
        reversed_path.push_back(current_state);
        State next_state;

        if (current_state == M) {
            int src = ws.get_ptrM(i, j);
            next_state = (src == 0) ? M : ((src == 1) ? I : D);
            i--; j--;
        } 
        else if (current_state == I) {
            int src = (j > 0) ? ws.get_ptrI(i, j) : ((i == 1) ? 0 : 1);
            next_state = (src == 0) ? M : I;
            i--;
        } 
        else {
            int src = ws.get_ptrD(i, j);
            next_state = (src == 0) ? M : D;
            j--;
        }
        current_state = next_state;
    }

    int start_j = j; 
    std::reverse(reversed_path.begin(), reversed_path.end());

    return {start_j, reversed_path};
}

void ProfileHMM::update_parameters_from_paths(const std::vector<std::string>& reads, 
                                              const std::vector<AlignmentPath>& paths) {
    float pseudo_count = 0.1f;
    std::vector<float> c_e_M((L + 1) * alpha_size, pseudo_count);
    std::vector<float> c_e_I((L + 1) * alpha_size, pseudo_count);

    for (size_t r = 0; r < reads.size(); ++r) {
        const std::string& read = reads[r];
        int i = 0, j = paths[r].start_j; 

        for (State curr_state : paths[r].states) {
            if (curr_state == M) {
                j++; i++;
                if ((size_t)i <= read.length() && j <= L) 
                    c_e_M[j * alpha_size + get_char_idx(read[i-1])] += 1.0f;
            } else if (curr_state == I) {
                i++;
                if ((size_t)i <= read.length() && j <= L) 
                    c_e_I[j * alpha_size + get_char_idx(read[i-1])] += 1.0f;
            } else if (curr_state == D) {
                j++;
            }
        }
    }

    for (int j = 1; j <= L; ++j) {
        float sum_e_M = 0, sum_e_I = 0;
        for (int k = 0; k < alpha_size; ++k) { 
            sum_e_M += c_e_M[j * alpha_size + k]; 
            sum_e_I += c_e_I[j * alpha_size + k]; 
        }

        for (int k = 0; k < alpha_size; ++k) {
            int idx = j * alpha_size + k;
            e_M[idx] = std::log(c_e_M[idx] / sum_e_M); 
            e_I[idx] = std::log(c_e_I[idx] / sum_e_I);
        }
    }
}

std::string ProfileHMM::extract_consensus(int min_depth, const std::vector<AlignmentPath>& final_paths) const {
    std::string consensus = "";
    std::vector<int> match_depth(L + 1, 0);
    
    for (const auto& path : final_paths) {
        int j = path.start_j; 
        for (State state : path.states) {
            if (state == M) { j++; if (j <= L) match_depth[j]++; }
            else if (state == D) { j++; }
        }
    }

    for (int j = 1; j <= L; ++j) {
        if (match_depth[j] < min_depth) continue; 

        float max_log_prob = LOG_ZERO;
        char best_char = alphabet[0];

        for (int k = 0; k < alpha_size; ++k) {
            if (e_M[j * alpha_size + k] > max_log_prob) {
                max_log_prob = e_M[j * alpha_size + k];
                best_char = alphabet[k];
            }
        }
        consensus += best_char;
    }
    return consensus;
}

std::string build_hmm_consensus(const std::string& backbone, 
                                const std::vector<std::string>& reads,
                                const std::vector<char>& valid_chars,
                                char ref_char,
                                int iterations,
                                int min_depth) {
    if (backbone.empty() || reads.empty() || valid_chars.empty()) return "";

    /*
    fprintf(stderr, "[telohmm] Backbone length: %zu; #reads: %zu; Min depth: %d\n", 
            backbone.length(), reads.size(), min_depth);
    */

    py::gil_scoped_release release;

    ProfileHMM hmm(backbone, valid_chars, ref_char);
    std::vector<AlignmentPath> all_paths(reads.size());
    
    ViterbiWorkspace ws;

    for (int iter = 0; iter < iterations; ++iter) {
        for (size_t r = 0; r < reads.size(); ++r) {
            all_paths[r] = hmm.viterbi_align(reads[r], ws);
        }
        hmm.update_parameters_from_paths(reads, all_paths);
    }

    std::string consensus = hmm.extract_consensus(min_depth, all_paths);

    return consensus;
}

PYBIND11_MODULE(telohmm, m) {
    m.doc() = "Optimized C++ Profile HMM for telomere TVR consensus";
    m.def("build_hmm_consensus", &build_hmm_consensus, 
          "Build consensus using Viterbi Training on Profile HMM",
          py::arg("backbone"), py::arg("reads"), py::arg("valid_chars"), py::arg("ref_char"),
          py::arg("iterations") = 3, py::arg("min_depth") = 3);
}
