#include <ctime>
#include "LSH.h"
#include "seqan/modifier.h"

LSH::LSH() :
    windowSize(100),
    windowOverlappingOffset(50),
    qGramSize(16),
    signatureLength(100),
    bandSize(2),
    limit(32),
    limit_skip(256),
    rngLimit(time(NULL)) {
}

LSH::LSH(unsigned long long seed) :
    LSH() {
    rngLimit.seed(seed);
}

void LSH::addChromosome(Chromosome &chromosome) {
    size_t end = seqan::length(chromosome.data) - windowSize;
    for (size_t i = 0; i < end; i += windowOverlappingOffset) {
        addReference(chromosome, i);
    }
}

LSH::ReferenceID LSH::addReference(Chromosome& chromosome, size_t beginPosition) {
    ReferenceID rid = references.size();
    references.push_back({&chromosome, beginPosition});
    auto it = seqan::begin(chromosome.data);
    std::unordered_set<unsigned> qgrams;
    createQGrams_iupac(it + beginPosition, it + (beginPosition + windowSize), qgrams);
    std::vector<unsigned> signatures;
    createSignatures(qgrams.begin(), qgrams.end(), signatures);
    addBandValues(signatures.begin(), signatures.end(), rid);
    return rid;
}

void LSH::findRead(ReadString read,
                   std::vector<ReferenceID>& out_referenceIDs,
                   std::vector<ReferenceID>& out_referenceIDs_revComplement) {
    findReadDirected(seqan::begin(read), seqan::end(read), out_referenceIDs);
    seqan::reverseComplement(read);
    findReadDirected(seqan::begin(read), seqan::end(read), out_referenceIDs_revComplement);
}


