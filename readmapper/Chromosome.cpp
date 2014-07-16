#include "Chromosome.h"

Chromosome::Chromosome(std::string name, ReferenceString data, std::vector<Variant> variants) :
    name(name),
    data(data),
    variants(variants) {
}

void Chromosome::addSNP(size_t pos, ReferenceChar snp) {
    assert(seqan::length(data) > pos);
    data[pos] = seqan::Iupac(snp).value | seqan::Iupac(data[pos]).value;
}

void Chromosome::addInsertion(size_t pos, ReferenceString insertionString) {
    assert(seqan::length(data) > pos);
    variants.emplace_back(pos, insertionString);
}

void Chromosome::addDeletion(size_t posAfterDeletion, size_t deletionLength) {
    assert(seqan::length(data) > posAfterDeletion);
    assert(posAfterDeletion >= deletionLength);
    variants.emplace_back(posAfterDeletion, deletionLength);
}
