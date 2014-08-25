#include "align/BacktracingMatrix.h"
#include "align/Variant.h"

BacktracingMatrix::BacktracingMatrix(int maxErrors) :
    matrix(std::vector<BacktracingColumn>(1, BacktracingColumn(0, maxErrors))),
    variantInColumn(1, Variant::none),
    positionInColumn(1, 0),
    offsetInColumn(1, 0),
    maxErrors(maxErrors){
    // initialize first column with scores (row with index 0 is init'ed with score 0 in column constructor)
    for (int i = 1; i <= maxErrors + 1; i++) {
        matrix[0].addScore(i);
    }
}

int BacktracingMatrix::getScore(size_t column, size_t row) const {
    return matrix[column].getScore(row);
}

AlignmentType BacktracingMatrix::getDir(size_t column, size_t row) const {
    return matrix[column].getDir(row);
}

int BacktracingMatrix::getBtColumn(size_t column, size_t row) const {
    return matrix[column].getBtColumn(row);
}

void BacktracingMatrix::addScore(size_t column, int score) {
    matrix[column].addScore(score);
}

void BacktracingMatrix::addBt(size_t column, int btColumn, AlignmentType dir) {
    matrix[column].addBt(btColumn, dir);
}

BacktracingColumn& BacktracingMatrix::getColumn(size_t column) {
    return matrix[column];
}

ptrdiff_t BacktracingMatrix::BacktracingMatrix::getHighestColumnIndex() const {
    return matrix.size() - 1;
}

const Variant& BacktracingMatrix::getVariantAt(size_t column) const {
    return variantInColumn[column];
}

int BacktracingMatrix::getPositionAt(size_t column) const {
    return positionInColumn[column];
}

int BacktracingMatrix::getOffsetAt(size_t column) const {
    return offsetInColumn[column];
}

void BacktracingMatrix::newColumn(size_t pos) {
    newColumn(pos, Variant::none, 0); // means "no variant"
}

void BacktracingMatrix::newColumn(size_t pos, const Variant &variant, int offset) {
    variantInColumn.push_back(variant);
    positionInColumn.push_back(pos);
    offsetInColumn.push_back(offset);
    matrix.push_back(BacktracingColumn(getHighestColumnIndex() + 1, maxErrors));
}
