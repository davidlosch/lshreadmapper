#ifndef PG583_REFERENCEREADER_H
#define PG583_REFERENCEREADER_H

#include "Chromosome.h"
#include "types.h"

#include <unordered_set>
#include <vector>
#include <string>

class ReferenceReader {
    std::unordered_map<std::string, size_t> chromosomeNameToIndex;
    std::unordered_set<std::string> chromosomeFilter;
    std::vector<Chromosome> &chromosomes;
public:
    ReferenceReader(std::vector<Chromosome> &chromosomes,
                    std::unordered_set<std::string> chromosomeFilter = std::unordered_set<std::string>());

    int readReferences(const std::vector<std::string> &referenceFiles);
private:
    int readReferenceFile(const std::string &referenceFile);
};

#endif // PG583_REFERENCEREADER_H
