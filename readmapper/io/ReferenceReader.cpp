#include "io/ReferenceReader.h"
#include "Chromosome.h"
#include "Logger.h"
#include "types.h"

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
        chromosomeNameToIndex[Chromosome::generateName(chromosomes[i].id)] = i;
    }
    for (auto &referenceFile : referenceFiles) {
        readReferenceFile(referenceFile);
    }
    auto &info = logger.info() << "FASTA: List of chromosomes in reference: ";
    for (auto &chromosome : chromosomes) {
        info << '"' << Chromosome::generateName(chromosome.id) << '"' << ' ';
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
        if (chromosomes.size() == chromosomeFilter.size()) {
            logger.info() << "FASTA: Chromosome filter set, all filtered chromosomes read"
                         ", omitting subsequent reference entries" << std::endl;
            break;
        }
        ++recordCounter;
        seqan::CharString id;
        seqan::CharString refGenome;
        int res = seqan::readRecord(id, refGenome, refStream);
        if (res != 0) {
            logger.err() << "FASTA: Could not read record number " << recordCounter <<
                            " in reference file: " << referenceFile << std::endl;
            return 1;
        }

        std::string name = Chromosome::generateName(seqan::toCString(id));
        if (name.size() == 0) {
            logger.err() << "FASTA: Could not deduce chromosome id for reference: \"" << id << "\"" << std::endl;
            continue;
        }
        auto chromosomeIndex = chromosomeNameToIndex.find(name);
        if (chromosomeIndex != chromosomeNameToIndex.end()) {
            auto warn = [&](const std::string &type, const std::string &value) {
                logger.warn() << "FASTA: Duplicate reference chromosome ignored. "
                              << type << ": \"" << value << "\"." << std::endl; };
            warn("Chromosome", name);
            warn("Previous id", chromosomes[chromosomeIndex->second].id);
            warn("Current id", seqan::toCString(id));
            continue;
        }
        if (!chromosomeFilter.empty() && chromosomeFilter.find(name) == chromosomeFilter.end()) {
            logger.info() << "FASTA: Chromosome filter set. "
                          << "Omitting chromosome: \"" << name << "\", id: \"" << id << "\"." << std::endl;
            continue;
        }
        chromosomeNameToIndex.emplace(name, chromosomes.size());
        chromosomes.emplace_back(seqan::toCString(id));
        seqan::assign(chromosomes.back().data, refGenome);
        logger.info() << "Read reference: " << id << std::endl;
    }
    return 0;
}
