#ifndef PG583_CHROMOSOME_H
#define PG583_CHROMOSOME_H

#include "align/Variant.h"
#include "types.h"

#include <string>

class Chromosome {
public:
    std::string name;
    ReferenceString data;
    std::vector<Variant> variants;

    Chromosome(std::string name, ReferenceString data = "", std::vector<Variant> variants = std::vector<Variant>());

    void addSNP(size_t pos, ReferenceChar snp);
    void addInsertion(size_t pos, ReferenceString insertionString);
    void addDeletion(size_t posAfterDeletion, size_t deletionLength);
};

#endif // PG583_CHROMOSOME_H
