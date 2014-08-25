#include "align/BacktracingColumn.h"

BacktracingColumn::BacktracingColumn() :
    column(1, DELETEDCHAR),
    score(1, 0),
    indexOfThisColumn(0),
    maxErrors(5){
}

BacktracingColumn::BacktracingColumn(size_t indexOfThisColumn, int maxErrors) :
    column(1, indexOfThisColumn - 1),   // first entry of this column (jump one column to the left
    recDir(1, DELETEDCHAR),               // and indicate, that a reference char has been skipped
    score(1, 0),
    indexOfThisColumn(indexOfThisColumn),
    maxErrors(maxErrors) {
}

int BacktracingColumn::getScore(size_t index) const {
    if (index > scoreSize()) {
        return maxErrors + 1;
    } else {
        return score[index];
    }
}

int BacktracingColumn::getBtColumn(size_t index) const {
    if (index > btSize()) {
        return indexOfThisColumn;
    } else {
        return column[index];
    }
}

AlignmentType BacktracingColumn::getDir(size_t index) const {
    if(index > btSize()) {
        return INSERTEDCHAR;
    } else {
        return recDir[index];
    }
}

void BacktracingColumn::addScore(int s) {
    score.push_back(s);
}

void BacktracingColumn::addBt(int btColumn, AlignmentType dir) {
    column.push_back(btColumn);
    recDir.push_back(dir);
}

size_t BacktracingColumn::scoreSize() const {
    return score.size() - 1;
}

size_t BacktracingColumn::btSize() const {
    return column.size() - 1;
}

