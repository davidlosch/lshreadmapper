#ifndef PG583_CHROMOSOME_H
#define PG583_CHROMOSOME_H

#include "align/Variant.h"
#include "types.h"

#include <string>
#include <vector>
#include <stddef.h>

class Chromosome {
public:
    std::string id;
    ReferenceString data;
    std::vector<Variant> variants;

    Chromosome(const std::string &id);
    Chromosome(const std::string &id, const ReferenceString &data);
    Chromosome(const std::string &id, const ReferenceString &data, const std::vector<Variant> &variants);

    void addSNP(size_t pos, ReferenceChar snp);
    void addInsertion(size_t pos, ReferenceString insertionString);
    void addDeletion(size_t posAfterDeletion, size_t deletionLength);

    static std::string generateName(const std::string &id);
};

#endif // PG583_CHROMOSOME_H
