#include <cmath>
#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>
#include "mex.h"

/*
 * Computes a weighted data set from common outcome permutations. 
 * Juri Kuronen (2018)
*/

class Key {
public:
    std::vector<uint32_t> comb;
    uint64_t primary;
    std::vector<uint64_t> overflow;
    Key(std::vector<uint32_t>&& c, uint64_t init) {
        comb = std::move(c);
        primary = init;
    }
    Key(std::vector<uint32_t>&& c, uint64_t init, std::vector<uint64_t>&& rest) {
        comb = std::move(c);
        primary = init;
        overflow = std::move(rest);
    }
    bool operator==(const Key& key2) const {
        if (overflow.size() != key2.overflow.size()) return false;
        if (primary == key2.primary) {
            for (uint32_t i = 0; i < overflow.size(); ++i) {
                if (overflow[i] != key2.overflow[i]) return false;
            }
            return true;
        }
        return false;
    }
};

Key get_key(const std::vector<uint32_t>& data_vec, const std::vector<uint32_t>& mb, uint32_t i, uint32_t d, uint32_t noc_max) {
    std::vector<uint32_t> comb(mb.size());
    for (uint32_t x = 0; x < mb.size(); ++x) comb[x] = data_vec[i * d + mb[x]];
    uint64_t primary = 1;
    uint32_t vals = 64 / noc_max - 1;
    for (uint32_t x = 0; x < std::min<uint32_t>(vals, mb.size()); ++x) {
        primary <<= noc_max;
        primary += comb[x];
    }
    if (mb.size() <= vals) return Key(std::move(comb), primary);
    std::vector<uint64_t> overflow;
    uint64_t key = (1 << noc_max) + data_vec[i * d + mb[vals]];
    for (uint32_t x = vals + 1; x < mb.size(); ++x) {
        if (x % vals == 0) {
            overflow.push_back(key);
            key = 1;
        }
        key <<= noc_max;
        key += comb[x];
    }
    overflow.push_back(key);
    return Key(std::move(comb), primary, std::move(overflow));
}

class Key_hash {
public:
    size_t operator()(const Key& key) const {
        std::hash<uint64_t> hasher;
        size_t hash = hasher(key.primary);
        for (uint64_t x : key.overflow) {
            hash ^= hasher(x) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        } 
        return hash; 
    }
};

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check arguments
    if (nlhs != 2 || nrhs != 4) mexErrMsgTxt("Wrong number of arguments.");
    if (!mxIsDouble(prhs[0])) mexErrMsgTxt("Argument 'data' must be double matrix.");
    if (!mxIsDouble(prhs[1]) && mxGetNumberOfElements(prhs[1]) != 1) mexErrMsgTxt("Argument 'node' must be scalar double.");
    if (!mxIsDouble(prhs[2]) && (mxGetM(prhs[2]) == 1 || mxGetN(prhs[2]) == 1)) mexErrMsgTxt("Argument 'mb' must be double vector.");
    if (!mxIsDouble(prhs[3]) && mxGetNumberOfElements(prhs[3]) != 1) mexErrMsgTxt("Argument 'noc_max' must be scalar double.");

    // Read variables
    uint32_t n = mxGetM(prhs[0]);
    uint32_t d = mxGetN(prhs[0]);
    uint32_t y = mxGetScalar(prhs[1]) - 1;
    uint32_t noc_max = 1 << ((uint32_t) std::ceil(std::log2(mxGetScalar(prhs[3]))));
    const double * const data_ptr = (double*) mxGetData(prhs[0]);
    const double * const mb_ptr = (double*) mxGetData(prhs[2]);
    std::vector<uint32_t> data(n * d);
    for (uint32_t i = 0; i < n * d; ++i) data[(i % n) * d + i / n] = data_ptr[i];
    std::vector<uint32_t> mb(mxGetNumberOfElements(prhs[2]));
    for (uint32_t i = 0; i < mb.size(); ++i) mb[i] = mb_ptr[i] - 1;

    // Calculate weights
    std::unordered_map<Key, uint32_t, Key_hash> weights;
    mb.push_back(y);
    for (uint32_t i = 0; i < n; ++i) {
        Key key = get_key(data, mb, i, d, noc_max);
        ++weights[key];
    }

    // Output weighted data set and weights
    plhs[0] = mxCreateDoubleMatrix(weights.size(), d, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(weights.size(), 1, mxREAL);
    double* out_ptr0 = mxGetPr(plhs[0]);
    double* out_ptr1 = mxGetPr(plhs[1]);
    uint32_t row = 0;
    for (const auto& c : weights) {
        out_ptr1[row] = c.second;
        const std::vector<uint32_t>& row_vec = c.first.comb;
        for (uint32_t i = 0; i < row_vec.size(); ++i) {
            out_ptr0[mb[i] * weights.size() + row] = row_vec[i];
        }
        ++row;
    }

}
