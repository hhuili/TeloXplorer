/*
   Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>
   Licensed under the GNU General Public License v3.0
*/

#ifndef TELOHMM_HPP
#define TELOHMM_HPP

#include <vector>
#include <string>
#include <cstdint>

enum State { M = 0, I = 1, D = 2 };

struct AlignmentPath {
    int start_j;
    std::vector<State> states;
};

struct ViterbiWorkspace {
    int max_N = 0;
    int max_L = 0;
    int L_plus_1 = 0;
    
    std::vector<float> V_M, V_I, V_D;
    std::vector<uint8_t> ptr_M, ptr_I, ptr_D; 

    void ensure_capacity(int N, int L);

    inline float& get_VM(int i, int j) { return V_M[i * L_plus_1 + j]; }
    inline float& get_VI(int i, int j) { return V_I[i * L_plus_1 + j]; }
    inline float& get_VD(int i, int j) { return V_D[i * L_plus_1 + j]; }
    
    inline uint8_t& get_ptrM(int i, int j) { return ptr_M[i * L_plus_1 + j]; }
    inline uint8_t& get_ptrI(int i, int j) { return ptr_I[i * L_plus_1 + j]; }
    inline uint8_t& get_ptrD(int i, int j) { return ptr_D[i * L_plus_1 + j]; }
};

class ProfileHMM {
private:
    char ref_char;
    int L;
    std::string backbone;
    std::vector<char> alphabet;
    int alpha_size;

    std::vector<int> char_to_idx_table; 

    std::vector<float> t_MM, t_MI, t_MD;
    std::vector<float> t_IM, t_II;
    std::vector<float> t_DM, t_DD;

    std::vector<float> e_M;
    std::vector<float> e_I;

    void allocate_matrices();
    void initialize_probabilities(const std::string& backbone);

    inline int get_char_idx(char c) const {
        return char_to_idx_table[static_cast<unsigned char>(c)];
    }

public:
    ProfileHMM(const std::string& backbone, const std::vector<char>& valid_chars, char ref_char);

    AlignmentPath viterbi_align(const std::string& read, ViterbiWorkspace& ws) const;

    void update_parameters_from_paths(const std::vector<std::string>& reads, 
                                      const std::vector<AlignmentPath>& paths);

    std::string extract_consensus(int min_depth, const std::vector<AlignmentPath>& final_paths) const;
};

std::string build_hmm_consensus(const std::string& backbone, 
                                const std::vector<std::string>& reads,
                                const std::vector<char>& valid_chars,
                                char ref_char,
                                int iterations,
                                int min_depth);

#endif
