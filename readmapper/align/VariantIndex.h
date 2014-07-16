#ifndef PG583_VARIANTINDEX_H
#define PG583_VARIANTINDEX_H

#include "align/Variant.h"

#include <vector>
#include <cstddef>
#include <assert.h>

struct Interval {
    size_t start;
    size_t end;
};

class VariantIndex {
public:
    VariantIndex(std::vector<Variant> &variantList, size_t highestPosition);
    unsigned int getVariantsAtPosition(size_t position) const;
    bool isVariantAtPosition(size_t position) const;
    const Variant& getVariant(size_t position) const;
    size_t getVariantCount() const;
private:
    std::vector<Variant> &variantList;
    std::vector<unsigned long long> existsVariant;
    std::vector<size_t> listPositions;
    std::vector<std::vector<size_t>> listOffset;
    unsigned int popcount(unsigned long long x) const;
};

#endif // PG583_VARIANTINDEX_H
