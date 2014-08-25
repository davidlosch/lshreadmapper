#ifndef PG583_IO_REFERENCE_READER_H
#define PG583_IO_REFERENCE_READER_H

#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <stddef.h>

class Chromosome;

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

#endif // PG583_IO_REFERENCE_READER_H
