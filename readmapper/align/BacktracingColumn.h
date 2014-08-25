#ifndef PG583_ALIGN_BACKTRACING_COLUMN_H
#define PG583_ALIGN_BACKTRACING_COLUMN_H

#include <vector>
#include <stddef.h>

enum AlignmentType : char {
    MATCH,          // read and reference chars match
    DELETEDCHAR,    // char from reference is not in read
    INSERTEDCHAR,   // char from read is not in reference
    MISMATCH,       // read and referece chars do not match
    NONE            // special case, used for deletion variants
};

class BacktracingColumn {
private:
    std::vector<int> column;            // column to jump to
    std::vector<AlignmentType> recDir;  // recursion direction
    std::vector<int> score;             // score at this field
    size_t indexOfThisColumn;
    int maxErrors;
public:
    BacktracingColumn();
    BacktracingColumn(size_t indexOfThisColumn, int maxErrors);
    int getScore(size_t index) const;
    int getBtColumn(size_t index) const;
    AlignmentType getDir(size_t index) const;
    void addScore(int s);
    void addBt(int btColumn, AlignmentType dir);
    size_t scoreSize() const;
    size_t btSize() const;
};

#endif // PG583_ALIGN_BACKTRACING_COLUMN_H
