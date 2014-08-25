#ifndef PG583_ALIGN_SEMI_GLOBAL_ALIGNER_H
#define PG583_ALIGN_SEMI_GLOBAL_ALIGNER_H

#include "align/VariantIndex.h"
#include "align/Alignment.h"
#include "Chromosome.h"
#include "types.h"

#include <seqan/sequence.h>

#include <vector>
#include <stddef.h>

class BacktracingMatrix;
class BacktracingColumn;

class SemiGlobalAligner {
private:
    /**
     * @brief variantIndex
     */
    VariantIndex variantIndex;
    typedef seqan::Infix<const decltype(((Chromosome*)nullptr)->data)>::Type ReferenceInfix;
    Alignment backtrace(ReferenceInfix &text,
                         const ReadString &pattern,
                         BacktracingMatrix &b,
                         int position);
    BacktracingColumn& alignInsertion(const Variant &variant,
                        const ReadString &pattern,
                        int maxError,
                        size_t lastRowArg,
                        size_t refPos,
                        int &bestMinimum,
                        std::vector<int> &bestMinimumPositions,
                        BacktracingMatrix &b);
    void printMatrix(const ReferenceInfix &text, const ReadString &pattern, const BacktracingMatrix &b) const;

public:
    SemiGlobalAligner(VariantIndex &variantIndex);
    std::vector<Alignment> align(const Chromosome &chromosome, size_t startPos, size_t endPos, const ReadString &pattern,
                    int maxError);
    int charDistance(seqan::Iupac p, seqan::Iupac t) const;
};

#endif // PG583_ALIGN_SEMI_GLOBAL_ALIGNER_H
