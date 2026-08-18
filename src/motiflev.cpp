/*
   Copyright (C) 2025 Huihui Li <hhui.li@outlook.com>
   Licensed under the GNU General Public License v3.0
*/

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using MotifToken = std::pair<std::string, int>;
using MotifSequence = std::vector<MotifToken>;

class Aligner {
private:
    double max_possible_cost(const MotifSequence& seq1, const MotifSequence& seq2) const {
        double cost = 0.0;
        for (const auto& t : seq1) cost += deletion_cost(t);
        for (const auto& t : seq2) cost += insertion_cost(t);
        return cost;
    }

public:
    double slippage_weight;  // Motif Match, Count differ (Replication Slippage)
    double mismatch_open;    // Motif Mismatch initial penalty (Substitution)
    double mismatch_ext;     // Motif Mismatch extension penalty (log-scaled)
    double gap_open;         // Gap Open Penalty
    double gap_extend;       // Gap Extension Penalty (log-scaled)

    Aligner(double slippage = 1.0, double mis_base = 3.0, double mis_ext = 0.5, 
            double g_open = 2.0, double g_extend = 0.5)
        : slippage_weight(slippage), mismatch_open(mis_base), mismatch_ext(mis_ext), 
          gap_open(g_open), gap_extend(g_extend) {}

    double insertion_cost(const MotifToken& b) const {
        return gap_open + gap_extend * std::log2(1.0 + b.second);
    }

    double deletion_cost(const MotifToken& a) const {
        return gap_open + gap_extend * std::log2(1.0 + a.second);
    }

    double alignment_cost(const MotifToken& a, const MotifToken& b) const {
        if (a.first == b.first) {
            return slippage_weight * std::log2(1.0 + std::abs(a.second - b.second));
        } else {
            return mismatch_open + mismatch_ext * std::log2(1.0 + a.second + b.second);
        }
    }

    double distance(const MotifSequence& seq1, const MotifSequence& seq2) const {
        if (seq1.empty() && seq2.empty()) return 0.0;

        int m = seq1.size();
        int n = seq2.size();

        std::vector<double> dp(n + 1, 0.0);

        for (int j = 1; j <= n; ++j) {
            dp[j] = dp[j-1] + insertion_cost(seq2[j-1]);
        }

        double current_del_cost = 0.0; 

        for (int i = 1; i <= m; ++i) {
            double prev_diag = dp[0]; 
            current_del_cost += deletion_cost(seq1[i-1]);
            dp[0] = current_del_cost;

            for (int j = 1; j <= n; ++j) {
                double temp = dp[j]; 
                double cost_del = dp[j] + deletion_cost(seq1[i-1]);
                double cost_ins = dp[j-1] + insertion_cost(seq2[j-1]);
                double cost_sub = prev_diag + alignment_cost(seq1[i-1], seq2[j-1]);

                dp[j] = std::min({cost_del, cost_ins, cost_sub});
                prev_diag = temp;
            }
        }

        double max_cost = max_possible_cost(seq1, seq2);
        if (max_cost <= 0.0) return 0.0;

        double norm_dist = dp[n] / max_cost;
        return std::min(1.0, std::max(0.0, norm_dist));
    }
};

namespace py = pybind11;

double default_distance(const MotifSequence& seq1, const MotifSequence& seq2) {
    Aligner default_aligner;
    return default_aligner.distance(seq1, seq2);
}

PYBIND11_MODULE(motiflev, m) {
    m.doc() = "Motif-aware Levenshtein distance calculation for tandem repeats and token sequences.";

    py::class_<Aligner>(m, "Aligner")
        .def(py::init<double, double, double, double, double>(),
             py::arg("slippage_weight") = 1.0,
             py::arg("mismatch_open") = 3.0,
             py::arg("mismatch_ext") = 0.5,
             py::arg("gap_open") = 2.0,
             py::arg("gap_extend") = 0.5,
             "Initialize the Motif Aligner with custom penalty parameters (e.g., replication slippage weight).")
        .def("distance", &Aligner::distance,
             py::arg("seq1"), py::arg("seq2"),
             "Calculate the normalized Motif-aware Levenshtein Distance [0.0, 1.0].");

    m.def("distance", &default_distance, 
          py::arg("seq1"), py::arg("seq2"),
          "Calculate the normalized Motif-aware Levenshtein distance using default parameters.");
}
