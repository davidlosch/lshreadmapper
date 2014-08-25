#include "align/SemiGlobalAligner.h"
#include "align/BacktracingMatrix.h"
#include "align/Variant.h"

#include <iostream>
#include <sstream>
#include <algorithm>

using std::vector;

SemiGlobalAligner::SemiGlobalAligner(VariantIndex &variantIndex) :
    variantIndex(variantIndex) {
}

int SemiGlobalAligner::charDistance(seqan::Iupac p, seqan::Iupac t) const {
    return 0 == (p.value & t.value);
}

vector<Alignment> SemiGlobalAligner::align(const Chromosome &chromosome, size_t startPos, size_t endPos,
                                   const ReadString &pattern, int maxError) {
    ReferenceInfix text(chromosome.data, startPos, endPos);
    // =========================================
    // 1. Calculate Matrix by Dynamic Programming
    // =========================================
    BacktracingMatrix b(maxError); // backtracing matrix
    size_t lastRow = std::min((size_t)maxError, seqan::length(pattern));
    int bestMinimum = maxError + 1;
    vector<int> bestMinimumPositions;

    for (size_t refPos = 0; refPos < seqan::length(text); refPos++) {
        // iterate over columns

        // optimization
        lastRow = std::min(lastRow, seqan::length(pattern));
        size_t rowCancel = lastRow;
        lastRow = 1;

        // === detect variants ===
        // vector to store relevant columns from all variants - including the last column of the current matrix
        vector<BacktracingColumn> columns(1, b.getColumn(b.getHighestColumnIndex()));
        vector<int> lastColumnOfVariant(1, b.getHighestColumnIndex()); // stored end position of each variant

        if (variantIndex.isVariantAtPosition(startPos + refPos)) {
            // if variants exist, iterate over all variants at current position
            unsigned int variantListStart = variantIndex.getVariantsAtPosition(startPos + refPos);
            unsigned int k = variantListStart;
            while (k < variantIndex.getVariantCount() && variantIndex.getVariant(k).getPosition() == startPos + refPos) {
                const Variant &currentVar = variantIndex.getVariant(k); // better use pointer?
                BacktracingColumn newColumn;
                switch(currentVar.getType()) {
                case(Variant::INSERTION):
                    newColumn = alignInsertion(currentVar, pattern, maxError, rowCancel, refPos, bestMinimum,
                                               bestMinimumPositions, b);
                    break;
                case(Variant::DELETION):
                    if (refPos - currentVar.getDeletionLength() > 0) {
                        b.newColumn(refPos, currentVar, 0);
                        for (size_t j = 1; j <= rowCancel; j++) {
                            b.addBt(b.getHighestColumnIndex(), b.getHighestColumnIndex()- 1 - currentVar.getDeletionLength(), NONE);
                        }
                        newColumn = b.getColumn(b.getHighestColumnIndex() - 1 - currentVar.getDeletionLength());
                    }
                    break;
                case(Variant::COMPLEX):
                    // not supported yet
                    break;
                }
                // this step is necessary, because columns from variants could be "longer" then the reference one
                rowCancel = std::max(rowCancel, newColumn.scoreSize());
                // insert new column into vector
                columns.push_back(newColumn);
                lastColumnOfVariant.push_back(b.getHighestColumnIndex());
                k++;
            }
        }

        // add new column to b
        b.newColumn(refPos);

        // fill score and backtracing column
        for (size_t row = 1; row <= rowCancel; row++) {
            // initialize with "downwards recursion" (= insertion)
            int recursionColumn = b.getHighestColumnIndex();
            int minimum = b.getScore(recursionColumn, row-1) + 1;
            AlignmentType recursionDirection = INSERTEDCHAR;

            // iterate over rows and check "matching" and "sidewards" recursion
            for (size_t k = 0; k < columns.size(); k++) {
                seqan::Iupac p = pattern[row-1];
                seqan::Iupac t = text[refPos];
                bool matched = charDistance(p, t) == 0;
                int toCheck = columns[k].getScore(row-1) + charDistance(p, t);
                if (toCheck < minimum) {
                    minimum = toCheck;
                    recursionDirection = matched ? MATCH : MISMATCH; // match
                    recursionColumn = lastColumnOfVariant[k];
                }
                toCheck = columns[k].getScore(row) + 1;
                if (toCheck < minimum) {
                    minimum = toCheck;
                    recursionDirection = DELETEDCHAR; // deletion
                    recursionColumn = lastColumnOfVariant[k];
                }
            }

            // fill columns
            b.addScore(b.getHighestColumnIndex(), minimum);
            b.addBt(b.getHighestColumnIndex(), recursionColumn, recursionDirection);

            // increment steps for next column
            if (minimum <= maxError)
                lastRow = row+1;

            // add/clear best positions
            if (row == seqan::length(pattern) && minimum <= bestMinimum) {
                if (row == seqan::length(pattern) && minimum < bestMinimum) {
                    bestMinimumPositions.clear();
                    bestMinimum = minimum;
                }
                bestMinimumPositions.push_back(b.getHighestColumnIndex());
            }
        }
    }

    // =========================================
    // 2. Calculate Alignment with Matrix
    // =========================================

//    printMatrix(text, pattern, b);

    vector<Alignment> alignments;
    if (bestMinimum <= maxError) {
        for (int position : bestMinimumPositions) {
            alignments.push_back(backtrace(text, pattern, b, position));
        }
    }
    return alignments;
}

Alignment SemiGlobalAligner::backtrace(ReferenceInfix &text, const ReadString &pattern, BacktracingMatrix &b, int position) {
    int currentColumn = position;
    int currentRow = seqan::length(pattern);
    Alignment alignment;
    while (currentRow > 0) {
        alignment.addUsedVariant(b.getVariantAt(currentColumn));
        AlignmentType recursionDirection = b.getDir(currentColumn, currentRow);
        int rowDecrease = 0;
        switch (recursionDirection) {
        case DELETEDCHAR:
            alignment.addCigarChar('D');
            if(&b.getVariantAt(currentColumn) == &Variant::none) {
                alignment.addAlignmentChar(text[b.getPositionAt(currentColumn)]);
            } else {
                alignment.addAlignmentChar(b.getVariantAt(currentColumn).getInsertionString()[b.getOffsetAt(currentColumn)]);
            }
            alignment.incrementErrorCount();
            break;
        case INSERTEDCHAR:
            alignment.addCigarChar('I');
            alignment.incrementErrorCount();
            rowDecrease++;
            break;
        case MATCH:
            alignment.addCigarChar('M');
            if(&b.getVariantAt(currentColumn) == &Variant::none) {
                alignment.addAlignmentChar(text[b.getPositionAt(currentColumn)]);
            } else {
                alignment.addAlignmentChar(b.getVariantAt(currentColumn).getInsertionString()[b.getOffsetAt(currentColumn)]);
            }
            rowDecrease++;
            break;
        case MISMATCH:
            alignment.addCigarChar('X');
            if(&b.getVariantAt(currentColumn) == &Variant::none) {
                alignment.addAlignmentChar(text[b.getPositionAt(currentColumn)]);
            } else {
                alignment.addAlignmentChar(b.getVariantAt(currentColumn).getInsertionString()[b.getOffsetAt(currentColumn)]);
            }
            alignment.incrementErrorCount();
            rowDecrease++;
            break;
        case NONE:
            break;
        default:
            break;
        }
        currentColumn = b.getBtColumn(currentColumn, currentRow);
        currentRow -= rowDecrease;
    }
    alignment.setPosition(b.getPositionAt(currentColumn));
    alignment.setStartVariant(b.getVariantAt(currentColumn));
    alignment.setStartOffset(b.getOffsetAt(currentColumn));
    return alignment;
}

BacktracingColumn& SemiGlobalAligner::alignInsertion(const Variant &variant,
                                       const ReadString &pattern,
                                       int maxError,
                                       size_t lastRowArg,
                                       size_t refPos,
                                       int &bestMinimum,
                                       vector<int> &bestMinimumPositions,
                                       BacktracingMatrix &b) {
    const auto &text = variant.getInsertionString();
    size_t lastRow = lastRowArg;

    for (size_t indelPos = 0; indelPos < seqan::length(text); indelPos++) {
        // iterate over columns

        // optimization
        size_t rowCancel = lastRow;
        lastRow = 1;

        int lastRowId = b.getHighestColumnIndex();
        b.newColumn(refPos, variant, indelPos);
        size_t currentColumn = b.getHighestColumnIndex();

        for (size_t row = 1; row <= rowCancel; row++) {
            // iterate over rows
            seqan::Iupac p = pattern[row-1];
            seqan::Iupac t = text[indelPos];
            bool matched = charDistance(p, t) == 0;
            AlignmentType recursionDirection = matched ? MATCH : MISMATCH; // match
            int minimum = b.getScore(currentColumn - 1, row - 1) + charDistance(p, t);
            if (b.getScore(currentColumn - 1, row) + 1 < minimum) {
                minimum = b.getScore(currentColumn - 1, row) + 1;
                recursionDirection = DELETEDCHAR; // deletion
            }
            if (b.getScore(currentColumn, row - 1) + 1 < minimum) {
                minimum = b.getScore(currentColumn, row -1 ) + 1;
                recursionDirection = INSERTEDCHAR; // insertion
            }

            // fill matrix
            b.addScore(currentColumn, minimum);
            b.addBt(currentColumn, lastRowId, recursionDirection);

            // increment steps for next column
            if (minimum <= maxError)
                lastRow = row+1;

            // add/clear best positions
            if (row == seqan::length(pattern) && minimum <= bestMinimum) {
                if (row == seqan::length(pattern) && minimum < bestMinimum) {
                    bestMinimumPositions.clear();
                    bestMinimum = minimum;
                }
                bestMinimumPositions.push_back(b.getHighestColumnIndex());
            }
        }
    }
    return b.getColumn(b.getHighestColumnIndex());
}

void SemiGlobalAligner::printMatrix(const ReferenceInfix &text, const ReadString &pattern, const BacktracingMatrix &b) const {
    const auto patternLength = seqan::length(pattern);
    const int maxIndex = b.getHighestColumnIndex();
    auto wout = [=]() -> std::ostream& { return std::cout << std::right << std::setw(2 + std::max(1, maxIndex / 10)); };

    wout() << ' ' << "=== Scores (Columns) ===" << std::endl;
    wout() << ' ';
    wout() << ' ';
    //TODO: wout() macht bei Texten irgendwie Probleme
    for (int j = 1; j <= b.getHighestColumnIndex(); j++) {
        if (&b.getVariantAt(j) == &Variant::none) {
            wout() << text[b.getPositionAt(j)];
        } else if (b.getVariantAt(j).getType() == Variant::INSERTION) {
            std::stringstream insertion;
            insertion << '(' << b.getVariantAt(j).getInsertionString()[b.getOffsetAt(j)] << ')';
            wout() << insertion.str();
        } else {
            wout() << '.';
        }
    }
    std::cout << std::endl;

    for (size_t i = 0; i <= patternLength; i++) {
        if (i > 0) {
            wout() << pattern[i-1];
        } else {
            wout() << ' ';
        }
        for (int j = 0; j <= maxIndex; j++) {
            wout() << b.getScore(j,i);
        }
        std::cout << std::endl;
    }

    wout() << ' ' << "=== Backtracing (Columns) ===" << std::endl;
    wout() << ' ';
    wout() << ' ';
    for (int j = 1; j <= maxIndex; j++) {
        wout() << j;
    }
    std::cout << std::endl;

    for (size_t i = 0; i <= patternLength; i++) {
        wout() << i;
        for (int j = 0; j <= maxIndex; j++) {
            wout() << b.getBtColumn(j,i);
        }
        std::cout << std::endl;
    }

    wout() << ' ' << "=== Backtracing (Direction) ===" << std::endl;
    wout() << ' ';
    wout() << ' ';
    for (int j = 1; j <= maxIndex; j++) {
        wout() << j;
    }
    std::cout << std::endl;

    for (size_t i = 0; i <= patternLength; i++) {
        wout() << i;
        for (int j = 0; j <= maxIndex; j++) {
            wout() << b.getDir(j,i);
        }
        std::cout << std::endl;
    }
}
