#ifndef PG583_BACKTRACINGMATRIX_H
#define PG583_BACKTRACINGMATRIX_H

#include "align/BacktracingColumn.h"
#include "align/Variant.h"

#include <vector>
#include <cstddef>
#include <functional>

class BacktracingMatrix {
private:
    std::vector<BacktracingColumn> matrix;
    std::vector<std::reference_wrapper<const Variant>> variantInColumn;
    std::vector<int> positionInColumn;
    std::vector<int> offsetInColumn;
    int maxErrors;
public:
    BacktracingMatrix(int maxErrors);
    int getScore(size_t column, size_t row) const;
    AlignmentType getDir(size_t column, size_t row) const;
    int getBtColumn(size_t column, size_t row) const;
    void addScore(size_t column, int score);
    void addBt(size_t column, int btColumn, AlignmentType dir);
    BacktracingColumn& getColumn(size_t column);
    ptrdiff_t getHighestColumnIndex() const;
    const Variant &getVariantAt(size_t column) const;
    int getPositionAt(size_t column) const;
    int getOffsetAt(size_t column) const;
    void newColumn(size_t pos);
    void newColumn(size_t pos, const Variant &variant, int offset);
};

#endif // PG583_BACKTRACINGMATRIX_H
