#ifndef PG583_IO_VARIANTS_READER_H
#define PG583_IO_VARIANTS_READER_H

#include "Logger.h"

#include <seqan/sequence.h>
#include <seqan/vcf_io.h>

#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>
#include <iostream>
#include <stddef.h>

class Chromosome;

class VariantsReader {
    enum class AltType {
        SUPPORTED, MISSING_ALLELE, BREAKEND, STRUCTURAL, UNKNOWN, MULTI_NUCLEOTIDE_POLYMORPHISM,
        SUBSTITUTIONS_FOLLOWED_BY_INDEL, SAME_AS_REFERENCE
    };

    std::unordered_map<std::string, size_t> chromosomeNameToIndex;
    std::unordered_set<std::string> chromosomeFilter;
    std::vector<Chromosome> &chromosomes;
public:
    VariantsReader(std::vector<Chromosome> &chromosomes,
                   std::unordered_set<std::string> chromosomeFilter = std::unordered_set<std::string>());

    int readVariants(const std::vector<std::string> &variantFiles);
private:
    int readVariantsFile(const std::string &variantFile);
    int processRecord(seqan::VcfRecord &record, const seqan::VcfIOContext &context);
    AltType processVariant(Chromosome &chromosome, size_t beginPos, const seqan::CharString &ref,
                           seqan::CharString &alts);
    AltType convertAltToUpperCase(seqan::CharString &alt);

    template<class TReader>
    int readVariantsWithReader(TReader &reader) {
        seqan::VcfHeader header;
        seqan::VcfIOContext context(header.sequenceNames, header.sampleNames);

        int res = seqan::read(header, reader, context, seqan::Vcf());
        if (res != 0) {
            logger.err() << "Reading VCF header. Code: " << res << std::endl;
            return 1;
        }

        int entryCount = 0;
        while (!seqan::atEnd(reader)) {
            ++entryCount;
            seqan::VcfRecord record;
            int res = seqan::readRecord(record, reader, context, seqan::Vcf());
            if (res != 0 && !seqan::atEnd(reader)) {
                logger.err() << "Reading VCF record. Error Code: " << res << ", Entry: " << entryCount << std::endl;
                return 1;
            }

            processRecord(record, context);
        }
        return 0;
    }
};

#endif // PG583_IO_VARIANTS_READER_H
