/*
   Copyright (C) 2025 Zepu Miao
   Licensed under the GNU General Public License v3.0
*/

/* Last Modified: 2026-03-15 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace py = pybind11;
using namespace std;

struct WinDiv { 
    string ass; 
    string tchr; 
    int start; 
    int pend; 
    float value; 
    int nnode; 
};

struct MnodeRan { int posStart; int posEnd; };

[[nodiscard]] float wMean(const vector<WinDiv> &winVec, int s, int e, int wmet) {
    double sumLen = 0.0, sumVal = 0.0;
    for (int v = s; v <= e; ++v) {
        double w = (wmet == 1) ? (winVec[v].pend - winVec[v].start + 1) : winVec[v].nnode;
        sumLen += w; sumVal += winVec[v].value * w;
    }
    return (sumLen == 0.0) ? 0.0f : static_cast<float>(sumVal / sumLen);
}

void storeValue(const vector<WinDiv> &winVec, int s, int e, int wmet, vector<WinDiv> &outVec) {
    int tlNode = 0;
    float val = (s == e) ? winVec[s].value : wMean(winVec, s, e, wmet);
    for (int v = s; v <= e; ++v) tlNode += winVec[v].nnode;
    outVec.push_back({winVec[s].ass, winVec[s].tchr, winVec[s].start, winVec[e].pend, val, (s == e ? winVec[s].nnode : tlNode)});
}

void mergeNode(vector<WinDiv> &vec, int dx, float dy, float base, int plateau, int wmet, vector<WinDiv> &outVec, int &lpos, int &rpos) {
    int n = (int)vec.size(), oneCount = 0;
    vector<uint8_t> flag(n, 0), maxFlag(n, 0);
    vector<MnodeRan> mran;

    while (oneCount < n) {
        float maxVal = -1e9; int k = -1;
        for (int i = 0; i < n; ++i) 
            if (!maxFlag[i] && (k == -1 || vec[i].value > maxVal)) { k = i; maxVal = vec[i].value; }
        
        if (k == -1 || maxVal <= base) break;
        maxFlag[k] = 1; oneCount++;
        float upL = maxVal + dy, lowL = maxVal - dy;
        int xstart = vec[k].start - dx, xend = vec[k].pend + dx;
        int rbound = k, lbound = k;
        bool rflag = false, lflag = false;
        float rpVal = 0, lpVal = 0;

        // Right Scan
        for (int j = k + 1, rpre = k; j < n; ++j) {
            if (vec[j].start > xend) break;
            if (vec[j].value >= lowL && vec[j].value <= upL) {
                xend = vec[j].pend + dx;
                if (flag[j]) break;
                if (rpre > k) {
                    float tval = wMean(vec, rpre + 1, j, wmet);
                    if (std::abs(tval - rpVal) <= dy) { rbound = j; rflag = true; } else break;
                } else {
                     if (vec[j].start - vec[rpre].pend > plateau) {
                        if ((vec[k].pend - vec[k].start + 1) >= plateau || (vec[j].pend - vec[j].start + 1) >= plateau) break;
                        if (j - rpre > 1) {
                            int mid = (vec[j].start + vec[rpre].pend) / 2, rnum = 0, lnum = 0;
                            for (int u = rpre + 1; u < j; ++u) { if (vec[u].pend <= mid) lnum++; else if (vec[u].start >= mid) rnum++; }
                            if (lnum == 0 && rnum == 0) { rbound = j; rpVal = wMean(vec, rpre, j, wmet); rflag = true; } else break;
                        } else break;
                     } else {
                        float tval = wMean(vec, rpre, j, wmet);
                        int lenCheck = ((vec[k].pend - vec[k].start + 1) >= plateau) ? 1 : ((vec[j].pend - vec[j].start + 1) >= plateau ? 2 : 0);
                        float diff = std::abs((lenCheck == 1 ? vec[rpre].value : (lenCheck == 2 ? vec[j].value : tval)) - tval);
                        if (lenCheck && diff > dy) break;
                        rbound = j; rpVal = tval; rflag = true;
                     }
                }
                if (rflag) rpre = rbound;
            } else if (flag[j]) break;
        }
        for (int s = k + 1; s <= rbound; ++s) if (!flag[s]) { flag[s] = maxFlag[s] = 1; oneCount++; }

        // Left Scan
        if (rflag) lpVal = rpVal;
        for (int t = k - 1, lpre = lbound; t >= 0; --t) {
            if (vec[t].pend < xstart) break;
            if (vec[t].value >= lowL && vec[t].value <= upL) {
                xstart = vec[t].start - dx;
                if (flag[t]) break;
                if (lpre < k) {
                    float tval = wMean(vec, t, lpre - 1, wmet);
                    if (std::abs(tval - lpVal) <= dy) { lbound = t; lflag = true; } else break;
                } else {
                    if (rflag) {
                        float tval = wMean(vec, t, lpre - 1, wmet);
                        if (std::abs(tval - lpVal) <= dy) { lbound = t; lflag = true; } else break;
                    } else {
                        if (vec[lpre].start - vec[t].pend > plateau) {
                            if ((vec[k].pend - vec[k].start + 1) >= plateau || (vec[t].pend - vec[t].start + 1) >= plateau) break;
                            if (lpre - t > 1) {
                                int mid = (vec[lpre].start + vec[t].pend) / 2, rnum = 0, lnum = 0;
                                for (int u = t + 1; u < lpre; ++u) { if (vec[u].pend <= mid) lnum++; else if (vec[u].start >= mid) rnum++; }
                                if (lnum == 0 && rnum == 0) { lbound = t; lpVal = wMean(vec, t, lpre, wmet); lflag = true; } else break;
                            } else if (vec[lpre].start - vec[t].pend > plateau) break;
                            else { lbound = t; lpVal = wMean(vec, t, lpre, wmet); lflag = true; }
                        } else {
                            float tval = wMean(vec, t, lpre, wmet);
                            int lenCheck = ((vec[k].pend - vec[k].start + 1) >= plateau) ? 1 : ((vec[t].pend - vec[t].start + 1) >= plateau ? 2 : 0);
                            float diff = std::abs((lenCheck == 1 ? vec[lpre].value : (lenCheck == 2 ? vec[t].value : tval)) - tval);
                            if (lenCheck && diff > dy) break;
                            lbound = t; lpVal = tval; lflag = true;
                        }
                    }
                }
                if (lflag) lpre = lbound;
            } else if (flag[t]) break;
        }
        for (int u = lbound; u < k; ++u) if (!flag[u]) { flag[u] = maxFlag[u] = 1; oneCount++; }
        
        if ((vec[k].pend - vec[k].start + 1) >= plateau || rflag || lflag) { flag[k] = 1; mran.push_back({lbound, rbound}); }
    }

    if (mran.empty()) { lpos = 0; rpos = n - 1; return; }

    sort(mran.begin(), mran.end(), [](const auto &a, const auto &b) { return a.posStart != b.posStart ? a.posStart < b.posStart : a.posEnd < b.posEnd; });

    lpos = (mran.back().posEnd < n - 1) ? mran.back().posEnd + 1 : mran.back().posStart;
    rpos = (mran.back().posEnd < n - 1) ? n - 1 : mran.back().posEnd;
    bool outLast = (mran.back().posEnd < n - 1);

    int prePos = 0;
    for (size_t r = 0; r < mran.size(); ++r) {
        int lowStart = prePos; bool lowFlag = true;
        if (prePos < mran[r].posStart) {
            for (int z = prePos; z < mran[r].posStart; ++z) {
                if (vec[z].value > base) {
                    if (lowFlag && z > lowStart) storeValue(vec, lowStart, z - 1, wmet, outVec);
                    storeValue(vec, z, z, wmet, outVec); lowFlag = false;
                } else if (!lowFlag) { lowStart = z; lowFlag = true; }
            }
            if (lowFlag) storeValue(vec, lowStart, mran[r].posStart - 1, wmet, outVec);
        }
        if (r < mran.size() - 1 || outLast) storeValue(vec, mran[r].posStart, mran[r].posEnd, wmet, outVec);
        prePos = mran[r].posEnd + 1;
    }
}

vector<WinDiv> scanNode(const vector<WinDiv> &inData, int win, int xdis, float ydis, float base, int plateau, int wmet) {
    vector<WinDiv> outData;
    outData.reserve(inData.size());
    vector<WinDiv> vec;
    string preAss = "", preChr = "";
    int count = 0;
    
    auto flush = [&](bool partial) {
        if (vec.empty()) return;
        int l = 0, r = 0;
        mergeNode(vec, xdis, ydis, base, plateau, wmet, outData, l, r);
        if (!partial) {
            storeValue(vec, l, r, wmet, outData);
            vec.clear(); count = 0;
        } else {
            vector<WinDiv> next; next.reserve(r - l + 1);
            for (int x = l; x <= r; ++x) next.push_back(vec[x]);
            vec = std::move(next); count = 0;
        }
    };

    for (const auto& item : inData) {
        if (preAss != "" && (item.ass != preAss || item.tchr != preChr)) flush(false);
        else if (count == win) flush(true);
        vec.push_back(item);
        preChr = item.tchr; preAss = item.ass; count++;
    }
    flush(false);
    return outData;
}

vector<WinDiv> mergeTop(const vector<WinDiv> &inData, int plateau, float topth, int wmet) {
    vector<WinDiv> outData;
    outData.reserve(inData.size());
    vector<WinDiv> vec; 
    WinDiv fDiv = {};
    string preAss = "", preChr = "";
    int midSum = 0;
    bool fFlag = false;

    auto flush = [&]() {
        if (fFlag) {
            outData.push_back(fDiv);
            for (const auto &d : vec) outData.push_back(d);
            vec.clear();
        }
        fFlag = false; midSum = 0;
    };

    for (const auto& item : inData) {
        int len = item.pend - item.start + 1;
        if (preAss != "" && (item.ass != preAss || item.tchr != preChr)) flush();
        
        if (!fFlag) {
            if (item.value < topth) outData.push_back(item);
            else { fDiv = item; fFlag = true; midSum = 0; }
        } else {
            if (item.value < topth) {
                midSum += len;
                if (midSum < plateau) vec.push_back(item);
                else {
                    outData.push_back(fDiv);
                    for (const auto &d : vec) outData.push_back(d);
                    outData.push_back(item);
                    fFlag = false; vec.clear();
                }
            } else {
                if (item.start - fDiv.pend <= plateau) {
                    double mVal = 0; int mLen = 0, mNode = 0;
                    for (const auto &d : vec) { mVal += d.value * (wmet==1? d.pend-d.start+1 : d.nnode); mLen += (d.pend-d.start+1); mNode += d.nnode; }
                    double fW = (wmet==1? fDiv.pend-fDiv.start+1 : fDiv.nnode), cW = (wmet==1? len : item.nnode);
                    fDiv.value = (float)((fDiv.value * fW + mVal + item.value * cW) / (fW + (wmet==1? mLen : mNode) + cW));
                    fDiv.pend = item.pend; fDiv.nnode += mNode + item.nnode;
                } else {
                    outData.push_back(fDiv);
                    for (const auto &d : vec) outData.push_back(d);
                    fDiv = item;
                }
                vec.clear(); midSum = 0;
            }
        }
        preChr = item.tchr; preAss = item.ass;
    }
    flush();
    return outData;
}

vector<WinDiv> finalNode(const vector<WinDiv> &inData, int xdis, float ydis, int plateau, int wmet) {
    vector<WinDiv> outData;
    outData.reserve(inData.size());
    vector<WinDiv> vec; 
    WinDiv fDiv = {};
    string preAss = "", preChr = "";
    int midSum = 0;
    bool fFlag = false;

    auto flush = [&]() {
        if (fFlag) {
            outData.push_back(fDiv);
            for (const auto &d : vec) outData.push_back(d);
            vec.clear();
        }
        fFlag = false; midSum = 0;
    };

    for (const auto& item : inData) {
        int len = item.pend - item.start + 1;
        if (preAss != "" && (item.ass != preAss || item.tchr != preChr)) flush();

        if (!fFlag) {
            if (len < plateau) outData.push_back(item);
            else { fDiv = item; fFlag = true; midSum = 0; }
        } else {
            if (len < plateau) {
                midSum += len;
                if (midSum < plateau) vec.push_back(item);
                else {
                    outData.push_back(fDiv);
                    for (const auto &d : vec) outData.push_back(d);
                    outData.push_back(item);
                    fFlag = false; vec.clear();
                }
            } else {
                bool merged = false;
                if (item.start - fDiv.pend <= xdis && std::abs(item.value - fDiv.value) <= ydis) {
                    double mVal = 0; int mLen = 0, mNode = 0;
                    for (const auto &d : vec) { mVal += d.value * (wmet==1? d.pend-d.start+1 : d.nnode); mLen += (d.pend-d.start+1); mNode += d.nnode; }
                    double fW = (wmet==1? fDiv.pend-fDiv.start+1 : fDiv.nnode), cW = (wmet==1? len : item.nnode);
                    double tv1 = (fDiv.value * fW + item.value * cW) / (fW + cW);
                    double tv2 = (fDiv.value * fW + mVal + item.value * cW) / (fW + (wmet==1? mLen : mNode) + cW);
                    
                    if (std::abs(tv1 - tv2) <= ydis) {
                        fDiv.pend = item.pend; fDiv.value = (float)tv2; fDiv.nnode += mNode + item.nnode; merged = true;
                    }
                }
                if (!merged) {
                    outData.push_back(fDiv);
                    for (const auto &d : vec) outData.push_back(d);
                    fDiv = item;
                }
                vec.clear(); midSum = 0;
            }
        }
        preChr = item.tchr; preAss = item.ass;
    }
    flush();
    return outData;
}

void getBase(const vector<WinDiv> &inData, float yper, float &baseVal, float &yDis, float &topth, float &maxAll) {
    if (inData.empty()) throw std::runtime_error("Input data is empty");
    vector<float> all;
    all.reserve(inData.size());
    for (const auto& d : inData) all.push_back(d.value);
    
    size_t n = all.size();
    auto getQ = [&](double p) { 
        auto it = all.begin() + (size_t)(n * p) - 1; 
        nth_element(all.begin(), it, all.end()); 
        return *it; 
    };
    
    baseVal = getQ(0.025);
    yDis = (getQ(0.975) - baseVal) * (yper / 100.0f);
    topth = getQ(0.95);
    maxAll = *max_element(all.begin(), all.end());
}

vector<WinDiv> run_bloom(
    const vector<WinDiv>& input_data,
    float ydis = -1.0f,
    float bvalue = -1.0f,
    float yper = 5.0f,
    int xdis = 500,
    int pwidth = 50,
    int iter = 5,
    int wmet = 1,
    bool mtop = false) 
{
    if (input_data.empty()) return {};

    float botVal=0, ryDis=0, topth=0, maxAll=0;
    getBase(input_data, yper, botVal, ryDis, topth, maxAll);

    float actual_bvalue = (bvalue < 0.0f) ? botVal : bvalue;
    float actual_ydis = (ydis < 0.0f) ? ryDis : ydis;

    vector<WinDiv> curData = input_data;
    
    if (mtop) {
        curData = mergeTop(curData, pwidth, topth, wmet);
    }

    size_t preCount = 0;
    for (int k = 0; k < iter; ++k) {
        curData = scanNode(curData, 100, xdis, actual_ydis, actual_bvalue, pwidth, wmet);
        if (k > 0 && curData.size() == preCount) break; // Early stopping
        preCount = curData.size();
    }

    return finalNode(curData, xdis, actual_ydis, pwidth, wmet);
}

PYBIND11_MODULE(bloom_core, m) {
    m.doc() = "BLOOM algorithm C++ extension for memory-efficient processing";

    py::class_<WinDiv>(m, "WinDiv")
        .def(py::init<std::string, std::string, int, int, float, int>())
        .def_readwrite("ass", &WinDiv::ass)
        .def_readwrite("tchr", &WinDiv::tchr)
        .def_readwrite("start", &WinDiv::start)
        .def_readwrite("pend", &WinDiv::pend)
        .def_readwrite("value", &WinDiv::value)
        .def_readwrite("nnode", &WinDiv::nnode)
        .def("__repr__", [](const WinDiv &w) {
            return "<WinDiv(ass='" + w.ass + 
                   "', tchr='" + w.tchr + 
                   "', start=" + std::to_string(w.start) + 
                   ", pend=" + std::to_string(w.pend) + 
                   ", value=" + std::to_string(w.value) + 
                   ", nnode=" + std::to_string(w.nnode) + ")>";
        });

    m.def("run_bloom", &run_bloom, "Run the BLOOM smoothing algorithm in memory",
          py::arg("input_data"),
          py::arg("ydis") = -1.0f,
          py::arg("bvalue") = -1.0f,
          py::arg("yper") = 5.0f,
          py::arg("xdis") = 500,
          py::arg("pwidth") = 50,
          py::arg("iter") = 5,
          py::arg("wmet") = 1,
          py::arg("mtop") = false,
          py::call_guard<py::gil_scoped_release>());
}
