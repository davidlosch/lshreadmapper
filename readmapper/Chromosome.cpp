#include "Chromosome.h"

#include <cassert>

Chromosome::Chromosome(const std::string &id) :
    Chromosome(id, "") {
}

Chromosome::Chromosome(const std::string &id, const ReferenceString &data) :
    Chromosome(id, data, std::vector<Variant>()) {
}

Chromosome::Chromosome(const std::string &id, const ReferenceString &data, const std::vector<Variant> &variants) :
    id(id),
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

std::string Chromosome::generateName(const std::string &id) {
    auto nameIt = seqan::begin(id);
    for (auto prefix : (const std::string[]) {"chrUn", "chr"}) {
        if (seqan::startsWith(id, prefix)) {
            nameIt += prefix.size();
            break;
        }
    }
    static const std::string ws = " \t";
    std::string name(nameIt, std::find_first_of(nameIt, seqan::end(id), ws.begin(), ws.end()));

    if (name == "M") { // Mitochondrium
        name = "MT";
    }

    return name;
}
