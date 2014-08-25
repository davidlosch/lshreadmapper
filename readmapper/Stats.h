#ifndef PG583_STATS_H
#define PG583_STATS_H

#include "module/VariantStatisticsModule.h"

#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <sstream>

void testStats();

template <class T>
std::string toString(T obj) {
     std::ostringstream ss;
     ss << obj;
     return ss.str();
}

// ======================================================================
// ====== WRITE
// ======================================================================
template <class Value>
void writeHistogramAsCSV(const std::vector<Value> &vector, std::ostream &strm, const std::string &headline = "") {
    std::cout << "writeHistogramAsCSV " << headline << std::endl;
    strm << headline << "\n";
    for (size_t i = 0; i < vector.size(); i++) {
        strm << i << ";" << vector[i] << "\n";
    }
    strm << std::endl;
}

template <class Index, class Value>
void writePairsAsCSV(const std::vector<std::pair<Index, Value>> &vector, std::ostream &strm,
                     const std::string &headline = "") {
    std::cout << "writePairsAsCSV " << headline << std::endl;
    strm << headline << "\n";
    for (size_t i = 0; i < vector.size(); i++) {
        strm << vector[i].first << ";" << vector[i].second << "\n";
    }
    strm << std::endl;
}

template <class Iterator>
void writePairsAsCSV(Iterator begin, Iterator end, std::ostream &strm, const std::string &headline = "") {
    std::cout << "writePairsAsCSV " << headline << std::endl;
    strm << headline << "\n";
    for (; begin != end; ++begin) {
        strm << begin->first << ";" << begin->second << "\n";
    }
    strm << std::endl;
}

template <class Iterable>
void writePairsAsCSV(const Iterable &it, std::ostream &strm, const std::string &headline = "") {
    writePairsAsCSV(begin(it), end(it), strm, headline);
}

template <class Value>
void writeMatrixAsCSV(const std::vector<std::vector<Value>> &vector, std::ostream &strm,
                      const std::string &headline = "") {
    std::cout << "writeMatrixAsCSV " << headline << std::endl;
    strm << headline << "\n";
    for (size_t i = 0; i < vector.size(); i++) {
        for (size_t j = 0; j < vector[i].size(); j++) {
            strm << vector[i][j] << ";";
        }
        strm << "\n";
    }
    strm << std::endl;
}

// ======================================================================
// ====== CONVERT
// ======================================================================
template <class Index, class Value, class ResultType = Value>
std::vector<std::pair<Index, ResultType>> convertMapToHistogram(const std::map<Index, Value> &map,
                                                                const std::vector<Index> &steps) {
    std::vector<std::pair<Index, ResultType>> result;
    result.reserve(steps.size() + 1);
    for (Index i: steps) {
        result.emplace_back(i, 0);
    }
    result.emplace_back(std::numeric_limits<Index>::max(), 0);

    size_t currentIdx = 0;
    for (std::pair<Index, Value> entry: map) {
        // map ist nach Index sortiert!
        while (currentIdx < result.size() && entry.first > result[currentIdx].first) {
            currentIdx++;
        }
        result[currentIdx].second += entry.second;
    }

    return result;
}

template <class ResultType, class Index, class Value>
std::vector<std::pair<Index, ResultType>> convertMapToHistogram(const std::map<Index, Value> &map,
                                                                const std::vector<Index> &steps, bool normalize) {
    std::vector<std::pair<Index, ResultType>> result = convertMapToHistogram<Index, Value, ResultType>(map, steps);
    if (normalize) {
        for (size_t i = 1; i < steps.size(); i++) {
            result[i].second /= steps[i] - steps[i-1];
        }
        result[0].second /= std::max(Index(1), steps.front() - map.begin()->first);
        result[steps.size()].second /= std::max(Index(1), map.rbegin()->first - steps.back());
    }
    return result;
}

// input and startPositions must be sorted!
template <class Pos, class Chromosome>
std::vector<std::pair<Chromosome, Pos>> convertPosToChromosomePos(const std::vector<Pos> &input,
                                            const std::vector<std::pair<Chromosome, Pos>> &startPositions) {
    std::vector<std::pair<Chromosome, Pos>> result(input.size());
    size_t curChromoIdx = 0;
    for (size_t i = 0; i < input.size(); i++) {
        size_t pos = input[i];
        auto nextChromoStartPos = [&](){return curChromoIdx+1 >= startPositions.size()
                                               ? std::numeric_limits<Pos>::max()
                                               : startPositions[curChromoIdx+1].second;};
        while (pos >= nextChromoStartPos()) {
            curChromoIdx++;
        }
        result[i].first = startPositions[curChromoIdx].first;
        result[i].second = pos - startPositions[curChromoIdx].second;
    }
    return result;
}

template <class Index, class Value>
std::vector<std::pair<Index, Value>> insertAllXvaluesToMap(const std::map<Index, Value> &map) {

    Index beginIndex = map.begin()->first;
    Index lng = map.rbegin()->first - map.begin()->first + 1;
    std::vector<std::pair<Index, Value>> result(lng);
    auto it = map.begin();
    for (Index i = 0; i < lng; i++) {
        Index idx = beginIndex + i;
        if (it->first == idx) {
            result[i] = std::pair<Index, Value>(idx, it->second);
            ++it;
        } else {
            result[i] = std::pair<Index, Value>(idx, 0);
        }
    }

    return result;
}


// ======================================================================
// ====== HISTOGRAM
// ======================================================================

template <class Pos, class Count>
std::map<size_t, size_t>  histogram_distanceBetweenVariants(const std::vector<std::pair<Pos, Count>> &sortedVector) {
    std::cout << "histogram_distanceBetweenVariants" << std::endl;
    std::map<size_t, size_t>  result;
    for (size_t i = 1; i < sortedVector.size(); i++) {
        size_t dist = sortedVector[i].first - sortedVector[i-1].first;
        if (result.find(dist) == result.end()) {
            result[dist] = 0;
        }
        result[dist]++;
    }
    return result;
}


template <class Pos, class Count>
std::vector<size_t> histogram_numberOfVariantsInQgram(const std::vector<std::pair<Pos, Count>> &sortedVector,
                                                      int q_length) {
    std::cout << "histogram_numberOfVariantsInQgram " << q_length << std::endl;
    std::vector<size_t> result(q_length+1, 0);
    size_t curBeginIdx = 0; // current begin index
    size_t curEndIdx = 0; // current end index
    size_t curPos = 0; // current position
    size_t curVariantNum = 0; // current number of variants
    while (curEndIdx < sortedVector.size()) {
        size_t nextBeginPos = curBeginIdx < sortedVector.size() ? sortedVector[curBeginIdx].first
                                                                : std::numeric_limits<size_t>::max();
        size_t nextEndPos = sortedVector[curEndIdx].first + q_length;
        if (nextBeginPos < nextEndPos) {
            size_t num = nextBeginPos - curPos;
            result[curVariantNum] += num;
            curPos = nextBeginPos;
            curVariantNum++;
            curBeginIdx++;
        } else {
            size_t num = nextEndPos - curPos;
            result[curVariantNum] += num;
            curPos = nextEndPos;
            curVariantNum--;
            curEndIdx++;
        }
    }
    return result;
}

template <class ResultType, class Pos, class Count>
std::vector<ResultType> numberOfCombinations(const std::vector<std::pair<Pos, Count>> &sortedVector, int q,
                                             const std::vector<ResultType> &limits) {
    std::vector<ResultType> result(limits.size(), 0);
    size_t curBeginIdx = 0; // current begin index
    size_t curEndIdx = 0; // current end index
    size_t curPos = 0; // current position
    size_t curVariantNum = 0; // current number of variants
    ResultType curCombinations = 1; // current number of Combinations
    while (curEndIdx < sortedVector.size()) {
        size_t nextBeginPos = curBeginIdx < sortedVector.size() ? sortedVector[curBeginIdx].first
                                                                : std::numeric_limits<size_t>::max();
        size_t nextEndPos = sortedVector[curEndIdx].first + q;
        if (nextBeginPos < nextEndPos) {
            size_t num = nextBeginPos - curPos;
            for (size_t i = 0; i < limits.size(); i++) {
                result[i] += num * std::min(curCombinations, limits[i]);
            }
            curPos = nextBeginPos;
            curVariantNum++;
            curCombinations *= (sortedVector[curBeginIdx].second + 1);
            curBeginIdx++;
        } else {
            size_t num = nextEndPos - curPos;
            for (size_t i = 0; i < limits.size(); i++) {
                result[i] += num * std::min(curCombinations, limits[i]);
            }
            curPos = nextEndPos;
            curVariantNum--;
            curCombinations /= (sortedVector[curEndIdx].second + 1);
            curEndIdx++;
        }
    }
    return result;
}

template <class ResultType, class Pos, class Count>
std::vector<ResultType> histogram_numberOfCombinations(const std::vector<std::pair<Pos, Count>> &sortedVector,
                                                       int q_max,
                                                       ResultType limit = std::numeric_limits<ResultType>::max()) {
    std::cout << "histogram_numberOfCombinations " << q_max << std::endl;
    std::vector<ResultType> limitVector{limit};
    std::vector<ResultType> result(q_max);
    result[0] = 0;
    for (size_t q = 1; q < result.size(); q++) {
        result[q] = numberOfCombinations(sortedVector, q, limitVector).front();
    }
    return result;
}

template <class ResultType, class Pos, class Count>
std::vector<std::vector<ResultType>> histogram_numberOfCombinations(
        const std::vector<std::pair<Pos, Count>> &sortedVector,
        const std::vector<int> &qValues,
        const std::vector<ResultType> &limitValues) {
    std::cout << "histogram_numberOfCombinations " << std::endl;
    std::vector<std::vector<ResultType>> result(qValues.size()+1);
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = std::vector<ResultType>(limitValues.size()+1);
        for (size_t j = 0; j < result[i].size(); j++) {
            result[i][j] = 0;
        }
    }
    for (size_t i = 0; i < qValues.size(); i++) {
        result[i+1][0] = qValues[i];
    }
    for (size_t j = 0; j < limitValues.size(); j++) {
        result[0][j+1] = limitValues[j];
    }

    for (size_t i = 0; i < qValues.size(); i++) {
        std::vector<ResultType> v = numberOfCombinations(sortedVector, qValues[i], limitValues);
        for (size_t j = 0; j < limitValues.size(); j++) {
            result[i+1][j+1] = v[j];
        }
    }
    return result;
}

template <class Pos, class Count>
std::map<Count, size_t> histogram_numberOfVariantsPerPosition(const std::vector<std::pair<Pos, Count>> &variants) {
    std::map<Count, size_t> result;
    for (auto entry: variants) {
        Count num = entry.second; // number of variants for this position
        if (result.find(num) == result.end()) {
            result[num] = 0;
        }
        result[num]++;
    }
    return result;
}

template <class Pos, class Count>
std::map<Count, size_t> histogram_indelLength(const std::vector<std::pair<Pos, Count>> &indelLength) {
    std::map<Count, size_t> result;
    for (auto entry: indelLength) {
        Count length = entry.second;
        if (result.find(length) == result.end()) {
            result[length] = 0;
        }
        result[length]++;
    }
    return result;
}

template <class Pos, class Count>
std::map<size_t, size_t> histogram_noGapRunLength(const std::vector<std::pair<Pos, Count>> &sortedVector,
                                                  Pos allowedGap = 0) {
    std::map<size_t, size_t> result;
    const size_t size = sortedVector.size();

    Pos lastPos = sortedVector[0].first;
    size_t curRunLength = 1;
    for (size_t i = 1; i <= size; i++) {
        Pos nextPos = i < size ? sortedVector[i].first : std::numeric_limits<Pos>::max();
        if (nextPos <= lastPos + allowedGap + 1) {
            curRunLength++;
        } else {
            if (result.find(curRunLength) == result.end()) {
                result[curRunLength] = 0;
            }
            result[curRunLength]++;
            curRunLength = 1;
        }
        lastPos = nextPos;
    }
    return result;
}

template <class ReturnType, class Pos, class Count>
std::vector<std::vector<ReturnType>> histogram_noGapRunLength(const std::vector<std::pair<Pos, Count>> &sortedVector,
                                                              const std::vector<Pos> &allowedGaps,
                                                              const std::vector<size_t> &steps) {
    std::vector<std::vector<ReturnType>> result(steps.size()+1);
    for (size_t i = 0; i < result.size(); i++) {
        result[i] = std::vector<ReturnType>(allowedGaps.size()+1);
        for (size_t j = 0; j < result[i].size(); j++) {
            result[i][j] = 0;
        }
    }
    for (size_t i = 0; i < steps.size(); i++) {
        result[i+1][0] = steps[i];
    }
    for (size_t j = 0; j < allowedGaps.size(); j++) {
        result[0][j+1] = allowedGaps[j];
    }

    for (size_t j = 0; j < allowedGaps.size(); j++) {
        std::vector<std::pair<Pos, ReturnType>> v =
                convertMapToHistogram<ReturnType>(histogram_noGapRunLength(sortedVector, allowedGaps[j]), steps, true);
        // the values larger than steps.back() are not copied!
        for (size_t i = 0; i < steps.size(); i++) {
            result[i+1][j+1] = v[i].second;
        }
    }
    return result;
}

template <class Pos, class Count>
Pos findOnePositionWithManyVariants(const std::vector<std::pair<Pos, Count>> &sortedVector, size_t windowSize,
                                    int variants) {
    size_t curBeginIdx = 0; // current begin index
    size_t curEndIdx = 0; // current end index
    size_t curPos = 0; // current position
    int curVariantNum = 0; // current number of variants
    while (curEndIdx < sortedVector.size()) {
        size_t nextBeginPos = curBeginIdx < sortedVector.size() ? sortedVector[curBeginIdx].first
                                                                : std::numeric_limits<size_t>::max();
        size_t nextEndPos = sortedVector[curEndIdx].first + windowSize;
        if (nextBeginPos < nextEndPos) {
            curPos = nextBeginPos;
            curVariantNum++;
            curBeginIdx++;
            if (curVariantNum >= variants) {
                // a position was found where enough variants are
                return curPos - windowSize + 1;
            }
        } else {
            curPos = nextEndPos;
            curVariantNum--;
            curEndIdx++;
        }

    }
    return std::numeric_limits<Pos>::max(); // not found
}

template <class Pos, class Count>
std::vector<Pos> findPositionsWithManyVariants(const std::vector<std::pair<Pos, Count>> &sortedVector, int windowSize,
                                               int variants, Pos jumpAfterFind = 0) {
    std::vector<Pos> result;
    size_t curBeginIdx = 0; // current begin index
    size_t curEndIdx = 0; // current end index
    size_t curPos = 0; // current position
    int curVariantNum = 0; // current number of variants
    while (curEndIdx < sortedVector.size()) {
        size_t nextBeginPos = curBeginIdx < sortedVector.size() ? sortedVector[curBeginIdx].first
                                                                : std::numeric_limits<size_t>::max();
        size_t nextEndPos = sortedVector[curEndIdx].first + windowSize;
        if (nextBeginPos < nextEndPos) {
            curPos = nextBeginPos;
            curVariantNum++;
            curBeginIdx++;
            if (curVariantNum >= variants) {
                // a position was found where enough variants are
                Pos p = curPos - windowSize + 1;
                if (result.size() == 0 || result.back() + jumpAfterFind < p) {
                    result.push_back(p);
                }
            }
        } else {
            curPos = nextEndPos;
            curVariantNum--;
            curEndIdx++;
        }

    }
    return result;
}

// return tuple of chromosome, position, length, number of variants
template <class Chromosome, class Pos, class Count>
std::vector<std::vector<std::string>> findPositionsWithManyVariants(
        const std::vector<std::pair<Pos, Count>> &sortedVector,
        const std::vector<std::pair<Chromosome, Pos>> &startPositions, size_t windowSize, int variants) {

    std::vector<std::vector<std::string>> result;
    enum {CHROMOSOME, POSITION, LENGTH, VARIANTS, ENUM_SIZE};
    size_t curBeginIdx = 0; // current begin index
    size_t curEndIdx = 0; // current end index
    size_t curPos = 0; // current position
    int curVariantNum = 0; // current number of variants
    size_t curChromoIdx = 0; // current index in startPositions
    long long lastStartPos = 123456;
    size_t lastBeginIdx = 123456;
    bool flag = false; // true if we are in a window with many variants
    auto nextChromoStartPos = [&](){return curChromoIdx+1 >= startPositions.size()
                                           ? std::numeric_limits<Pos>::max()
                                           : startPositions[curChromoIdx+1].second;};
    while (curEndIdx < sortedVector.size()) {
        size_t nextBeginPos = curBeginIdx < sortedVector.size() ? sortedVector[curBeginIdx].first
                                                                : std::numeric_limits<size_t>::max();
        size_t nextEndPos = sortedVector[curEndIdx].first + windowSize;
        if (nextBeginPos < nextEndPos) {
            curPos = nextBeginPos;
            curVariantNum++;
            curBeginIdx++;
            if (!flag && curVariantNum >= variants && nextBeginPos != nextEndPos) {
                // a position was found where enough variants are
                flag = true;
                lastStartPos = curPos - windowSize + 1;
                lastBeginIdx = curBeginIdx - curVariantNum;
                result.emplace_back(ENUM_SIZE);
                while (curPos >= nextChromoStartPos()) {
                    curChromoIdx++;
                }
                result.back()[CHROMOSOME] = toString(startPositions[curChromoIdx].first);
                result.back()[POSITION] = toString(lastStartPos - startPositions[curChromoIdx].second);
            }
        } else {
            curPos = nextEndPos;
            curVariantNum--;
            if (flag && curVariantNum < variants && nextBeginPos != nextEndPos) { // && lastStartPos + windowSize < curPos - windowSize) {
                flag = false;
                result.back()[LENGTH] = toString(curPos - lastStartPos);
                result.back()[VARIANTS] = toString(curBeginIdx - lastBeginIdx);
            }
            curEndIdx++;
        }

    }
    return result;
}

template <class Pos, class Count>
std::vector<std::pair<Count, Count>> manyVariantsExampleWindow(const std::vector<std::pair<Pos, Count>> &snps,
                                                               const std::vector<std::pair<Pos, Count>> &indels,
                                                               Pos position, size_t windowSize) {
    std::vector<std::pair<Count, Count>> result(windowSize, std::pair<Count, Count>(0,0));

    Pos endPos = position + windowSize;
    auto f = [](std::pair<Pos, Count> x, Pos p){return x.first < p;};
    size_t idxSNP = std::lower_bound(begin(snps), end(snps), position, f) - begin(snps);
    size_t idxIndels = std::lower_bound(begin(indels), end(indels), position, f) - begin(indels);
    while (idxSNP < snps.size() && snps[idxSNP].first < endPos) {
        result[snps[idxSNP].first - position].first += snps[idxSNP].second;
        idxSNP++;
    }
    while (idxIndels < indels.size() && indels[idxIndels].first < endPos) {
        result[indels[idxIndels].first - position].second += indels[idxIndels].second;
        idxIndels++;
    }

    return result;
}

//template <class Pos, class Count, class Length>
//void calculateStatistics(std::vector<std::pair<Pos, Count>> &snps,
//                         std::vector<std::pair<Pos, Count>> &indels,
//                         std::vector<std::pair<Pos, Length>> &indelLength) {

void calculateStatistics(VariantStatisticsModule &vsm);



#endif // PG583_STATS_H
