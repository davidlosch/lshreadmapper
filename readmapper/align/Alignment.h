#ifndef PG583_ALIGN_ALIGNMENT_H
#define PG583_ALIGN_ALIGNMENT_H

#include "types.h"

#include <string>
#include <vector>
#include <set>
#include <functional>
#include <stddef.h>

class Variant;

class Alignment {
public:
    Alignment();
    std::string getCigar();
    ReferenceString getAlignmentString();
    const std::set<std::reference_wrapper<const Variant>>& getUsedVariants();
    size_t getPosition();
    const Variant& getStartVariant();
    int getStartOffset();
    size_t getErrorCount();
    void addCigarChar(char c);
    void addAlignmentChar(ReferenceChar c);
    void addUsedVariant(const Variant &var);
    void setPosition(size_t pos);
    void setStartVariant(const Variant &var);
    void setStartOffset(int off);
    void incrementErrorCount();
private:
    // cigar string for alignment
    std::string cigar;
    // used alignment string
    ReferenceString alignmentString;
    // the list of variants used for this alignment
    std::set<std::reference_wrapper<const Variant>> usedVariants;
    // start index in reference
    size_t position;
    // the variant in which the alignments starts (or Variant::none if it starts in reference)
    std::reference_wrapper<const Variant> startVariant;
    // the position in this variant
    int startOffset;
    size_t errorCount;
    bool cigarInvalid;
    bool alignmentStringInvalid;
    std::string reversedCigar;
    ReferenceString reversedAlignmentString;
    int currentCharCounter;
    char currentChar;
    void flush();
};

#endif // PG583_ALIGN_ALIGNMENT_H
