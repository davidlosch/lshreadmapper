#include "io/ReferenceReader.h"
#include "Logger.h"

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

ReferenceReader::ReferenceReader(std::vector<Chromosome> &chromosomes,
                                 std::unordered_set<std::string> chromosomeFilter) :
    chromosomeFilter(chromosomeFilter),
    chromosomes(chromosomes) {
}

int ReferenceReader::readReferences(const std::vector<std::string> &referenceFiles) {
    logger.info() << "FASTA: Reading reference files." << std::endl;
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        chromosomeNameToIndex[chromosomes[i].name] = i;
    }
    for (auto &referenceFile : referenceFiles) {
        readReferenceFile(referenceFile);
    }
    auto &info = logger.info() << "FASTA: List of read chromosomes in reference: ";
    for (auto &chromosome : chromosomes) {
        info << '"' << chromosome.name << '"' << ' ';
    }
    info << std::endl;
    return 0;
}

int ReferenceReader::readReferenceFile(const std::string &referenceFile) {
    logger.info() << "Reading reference file: " << referenceFile << std::endl;
    seqan::SequenceStream refStream(referenceFile.c_str());
    if (!isGood(refStream)) {
        logger.err() << "FASTA: Could not open reference file: " << referenceFile << std::endl;
        return 1;
    }
    size_t recordCounter = 0;
    while (!seqan::atEnd(refStream)) {
        ++recordCounter;
        seqan::CharString id;
        seqan::CharString refGenome;
        int res = seqan::readRecord(id, refGenome, refStream);
        if (res != 0) {
            logger.err() << "FASTA: Could not read record " << recordCounter << " in reference file: " << referenceFile
                       << std::endl;
            return 1;
        }
        logger.info() << "Read reference: " << id << std::endl;
        static const seqan::CharString CHROMOSOME_TAG = "chr";
        if (!seqan::startsWith(id, CHROMOSOME_TAG)) {
            logger.err() << "FASTA: Could not deduce chromosome id for reference: " << id << std::endl;
            logger.err() << "FASTA: Chromosome id for CHROMOSOME_NAME must start with: " <<
                          '"' << CHROMOSOME_TAG << "CHROMOSOME_NAME " << '"' << std::endl;
            continue;
        }
        auto chromosomeNameItBegin = seqan::begin(id) + seqan::length(CHROMOSOME_TAG);
        std::string chromosomeName(chromosomeNameItBegin, std::find(chromosomeNameItBegin, seqan::end(id), ' '));
        auto chromosomeIndex = chromosomeNameToIndex.find(chromosomeName);
        if (chromosomeIndex != chromosomeNameToIndex.end()) {
            logger.err() << "FASTA: Duplicate reference chromosome ignored. Chromosome: " << chromosomeName << std::endl;
            continue;
        }
        if (!chromosomeFilter.empty()) {
            if (chromosomeFilter.find(chromosomeName) == chromosomeFilter.end()) {
                logger.info() << "FASTA: Chromosome filter set, ignored chromosome: " << chromosomeName << std::endl;
                continue;
            }
        }
        auto &chromosome = *chromosomes.emplace(chromosomes.end(), chromosomeName);
        seqan::assign(chromosome.data, refGenome);
        if (chromosomes.size() == chromosomeFilter.size()) {
            logger.info() << "FASTA: Chromosome filter set, all filtered chromosomes read"
                         ", omitting subsequent reference entries" << std::endl;
            return 0;
        }
    }
    return 0;
}
