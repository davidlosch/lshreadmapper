#include "io/VariantsReader.h"
#include "Chromosome.h"
#include "types.h"

#include <type_traits>

VariantsReader::VariantsReader(std::vector<Chromosome> &chromosomes, std::unordered_set<std::string> chromosomeFilter) :
    chromosomeFilter(chromosomeFilter),
    chromosomes(chromosomes) {
}

template<class TSeq1, class TSeq2, class TLen>
TLen commonPrefixLength(TSeq1 &&seq1, TSeq2 &&seq2, TLen minLen) {
    assert(minLen <= seqan::length(seq1));
    assert(minLen <= seqan::length(seq2));
    auto itBeginSeq1 = seqan::begin(seq1);
    auto dist = std::mismatch(itBeginSeq1, itBeginSeq1 + minLen, seqan::begin(seq2)).first - itBeginSeq1;
    assert(dist >= 0);
    assert((TLen) dist <= minLen);
    return (TLen) dist;
}

template<class TSeq1, class TSeq2>
std::pair<typename seqan::Infix<TSeq1>::Type, typename seqan::Infix<TSeq2>::Type>
trimCommonPrefixAndSuffix(TSeq1 &seq1, TSeq2 &seq2) {
    auto minLen = std::min(seqan::length(seq1), seqan::length(seq2));
    auto prefixLen = commonPrefixLength(seq1, seq2, minLen);
    auto suffixLen = (prefixLen == minLen) ? 0 :
                       commonPrefixLength(seqan::reverseString(seq1), seqan::reverseString(seq2), minLen - prefixLen);
    return std::make_pair(seqan::infix(seq1, prefixLen, seqan::length(seq1) - suffixLen),
                          seqan::infix(seq2, prefixLen, seqan::length(seq2) - suffixLen));
}

int VariantsReader::readVariants(const std::vector<std::string> &variantFiles) {
    logger.info() << "VCF: Reading variant files." << std::endl;
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        chromosomeNameToIndex[Chromosome::generateName(chromosomes[i].id)] = i;
    }
    for (auto &variantFile : variantFiles) {
        readVariantsFile(variantFile);
    }
    for (auto &chromosome : chromosomes) {
        std::sort(chromosome.variants.begin(), chromosome.variants.end());
    }
    return 0;
}

int VariantsReader::readVariantsFile(const std::string &variantFile) {
    const size_t bufferSize = 65536;
    // TODO: maybe check for magic number instead of file extension
    if (seqan::endsWith(variantFile, ".gz")) {
    #if SEQAN_HAS_ZLIB
        seqan::Stream<seqan::GZFile> f;
        if (!seqan::open(f, variantFile.c_str(), "rb")) {
            logger.err() << "VCF: GZip file has the wrong format! At file: " << variantFile << std::endl;
            return 1;
        }

        seqan::RecordReader<seqan::Stream<seqan::GZFile>, seqan::SinglePass<>> reader(f, bufferSize);
        readVariantsWithReader(reader);
    #else // SEQAN_HAS_ZLIB
        logger.err() << "VCF: ZLIB not available! At file: " << variantFile << std::endl;
        return 1;
    #endif // SEQAN_HAS_ZLIB
    } else {
        std::ifstream f(variantFile, std::ios::in | std::ios::binary);
        if (!f) {
            logger.err() << "VCF: Could not open \"" << variantFile << "\"." << std::endl;
            return 1;
        }

        seqan::RecordReader<std::ifstream, seqan::SinglePass<>> reader(f, bufferSize);
        readVariantsWithReader(reader);
    }
    return 0;
}

int VariantsReader::processRecord(seqan::VcfRecord &record, const seqan::VcfIOContext &context) {
    std::string chromosomeName = seqan::toCString((*context.sequenceNames)[record.rID]);

    if (!chromosomeFilter.empty() && chromosomeFilter.find(chromosomeName) == chromosomeFilter.end()) {
        return 0;
    }

    auto chromosomeIt = chromosomeNameToIndex.find(chromosomeName);
    if (chromosomeIt == chromosomeNameToIndex.end()) {
        logger.warn() << "VCF: Skipping entry. Chromosome not in reference: \"" << chromosomeName << '"' << std::endl;
        return 1;
    }
    auto &chromosome = chromosomes[chromosomeIt->second];

    seqan::toUpper(record.ref);
    seqan::StringSet<seqan::CharString> alts;
    seqan::splitString(alts, record.alt, ',');
    for (const auto &orgAlt : alts) {
        seqan::CharString alt = orgAlt;
        auto ret = processVariant(chromosome, record.beginPos, record.ref, alt);
        if (ret != AltType::SUPPORTED) {
            auto &log = logger.warn() << "VCF: ";
            switch(ret) {
            case AltType::SUPPORTED:
                assert(false); // not reachable
                break;
            case AltType::MISSING_ALLELE:
                log << "Missing allele '*' (due to upstream deletion) unsupported";
                break;
            case AltType::BREAKEND:
                log << "Breakend variants unsupported";
                break;
            case AltType::STRUCTURAL:
                log << "Structural variants / symbolic alleles unsupported";
                break;
            case AltType::UNKNOWN:
                log << "Unknown ALT-format";
                break;
            case AltType::SUBSTITUTIONS_FOLLOWED_BY_INDEL:
                log << "Substitution + indel variant not supported yet";
                break;
            case AltType::MULTI_NUCLEOTIDE_POLYMORPHISM:
                log << "Multi-nucleotide polymorphisms not supported yet";
                break;
            case AltType::SAME_AS_REFERENCE:
                log << "No variation, alt == ref";
                break;
            }
            log << ", ignored. At" <<
                   " Chromosome: \"" << chromosomeName << "\"" <<
                   " Position: " << record.beginPos <<
                   " Ref: \"" << record.ref << "\"" <<
                   " Alt: \"" << orgAlt << "\"" <<
            std::endl;
        }
    }
    return 0;
}
VariantsReader::AltType VariantsReader::processVariant(Chromosome &chromosome, size_t beginPos,
        const seqan::CharString &ref, seqan::CharString &alt) {
    // TODO: check if dbSNP RV-entries are not reversed:
    //       ##INFO=<ID=RV,Number=0,Type=Flag,Description="RS orientation is reversed">
    // TODO: handle telomeres (pos = 0, pos = N+1
    // TODO: maybe support: VCFv4.2.pdf: 'If any of the ALT alleles is a symbolic allele (an angle-bracketed
    //                                   ID String "<ID>") then the padding base is required and POS denotes the
    //                                   coordinate of the base preceding the polymorphism.'
    auto ret = convertAltToUpperCase(alt);
    if (AltType::SUPPORTED != ret) {
        return ret;
    }
    seqan::Infix<decltype(ref)>::Type infixRef;
    seqan::Infix<decltype(alt)>::Type infixAlt;
    std::tie(infixRef, infixAlt) = trimCommonPrefixAndSuffix(ref, alt);
    auto lengthRef = seqan::length(infixRef);
    auto lengthAlt = seqan::length(infixAlt);
    auto pos = beginPos + seqan::beginPosition(infixRef);
    if (lengthRef == lengthAlt) {
        if (lengthRef == 0) {
            // should not happen : alt == ref
            return AltType::SAME_AS_REFERENCE;
        } else if (lengthRef == 1) {
            chromosome.addSNP(pos, infixAlt[0]);
            return AltType::SUPPORTED;
        } else {
            return AltType::MULTI_NUCLEOTIDE_POLYMORPHISM;
        }
    } else if (lengthRef == 0) {
        // plain insertion
        ReferenceString insertionString;
        seqan::assign(insertionString, infixAlt);
        chromosome.addInsertion(pos, insertionString);
        return AltType::SUPPORTED;
    } else if (lengthAlt == 0) {
        // plain deletion
        chromosome.addDeletion(pos + seqan::length(infixRef), seqan::length(infixRef));
        return AltType::SUPPORTED;
    } else {
        return AltType::SUBSTITUTIONS_FOLLOWED_BY_INDEL;
    }
}

VariantsReader::AltType VariantsReader::convertAltToUpperCase(seqan::CharString &alt) {
    for (auto &c : alt) {
        switch (c) {
        case 'a':
        case 'c':
        case 'g':
        case 't':
        case 'n':
            c = std::toupper(c);
        case 'A':
        case 'C':
        case 'G':
        case 'T':
        case 'N':
            continue;
        case '*':
            return AltType::MISSING_ALLELE;
        case '[':
        case ']':
            return AltType::BREAKEND;
        case '>':
        case '<':
            return AltType::STRUCTURAL;
        default:
            return AltType::UNKNOWN;
        }
    }
    return AltType::SUPPORTED;
}
