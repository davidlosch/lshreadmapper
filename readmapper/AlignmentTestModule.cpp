#include "AlignmentTestModule.h"
#include "align/SemiGlobalAligner.h"
#include "io/ReferenceReader.h"
#include "io/VariantsReader.h"
#include "Chromosome.h"

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <tuple>


AlignmentTestModule::AlignmentTestModule() {

}

void AlignmentTestModule::addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName) {
    using namespace seqan;

    addUsageLine(argParser, moduleName + " -V <VCF_FILE> -r <REFERENCE_FASTA_FILES> -m <CHROMOSOME,START_POS,READ>");

    addSection(argParser, moduleName + " Options");

    ArgParseOption variantsFileOpt("V", "variants", "VCF file to use.", ArgParseOption::INPUTFILE, "variants.vcf",
                                   true);
    addOption(argParser, variantsFileOpt);

    ArgParseOption referenceFilesOpt("r", "reference", "Fasta files to use as reference", ArgParseOption::INPUTFILE,
                                     "reference.fa.gz", true);
    addOption(argParser, referenceFilesOpt);

    ArgParseOption mappedReadEntriesOpt("m", "mapped-reads", "Fasta files to use as reference", ArgParseOption::STRING,
                                  "ChromosomeID1,StartPosInChromosome1,ReadString1 "
                                  "ChromosomeID2,StartPosInChromosome2,ReadString2 ...", true);
    addOption(argParser, mappedReadEntriesOpt);

    ArgParseOption maxErrorsOpt("e", "max-errors", "Maximum number of tolerated errors", ArgParseOption::INTEGER, "N");
    setDefaultValue(maxErrorsOpt, 3);
    addOption(argParser, maxErrorsOpt);

    ArgParseOption refChromosomeFilterOpt("f", "filter", "Use only given chromosomes", ArgParseOption::STRING, "N",
                                          true);
    addOption(argParser, refChromosomeFilterOpt);
}


void AlignmentTestModule::readOptions(seqan::ArgumentParser &argParser) {
    seqan::getOptionValue(maxErrors, argParser, "e");
    auto chromosomes = seqan::getOptionValues(argParser, "f");
    chromosomeFilter.insert(chromosomes.begin(), chromosomes.end());
    variantFiles = seqan::getOptionValues(argParser, "V");
    referenceFiles = seqan::getOptionValues(argParser, "r");
    std::vector<std::string> mappedReadEntries(seqan::getOptionValues(argParser, "m"));
    for(std::string &readEntry : mappedReadEntries) {
        std::replace(readEntry.begin(), readEntry.end(), ',', ' ');
        std::istringstream ss(readEntry);
        std::string chromosome;
        if (!(ss >> chromosome)) {
            logger.err() << "Parsing read entry (chromosome): " << readEntry << std::endl;
            std::exit(1); // TODO: ...
        }
        size_t startPos;
        if (!(ss >> startPos)) {
            logger.err() << "Parsing read entry (start position): " << readEntry << std::endl;
            std::exit(1); // TODO: ...
        }
        std::string read;
        if (!(ss >> read)) {
            logger.err() << "Parsing read entry (read): " << readEntry << std::endl;
            std::exit(1); // TODO: ...
        }
    }
}

int AlignmentTestModule::run() {
    std::vector<Chromosome> chromosomes;

    Chromosome chromosome("testChromosome", "AAAGGGACCAAA", std::vector<Variant>());
    chromosome.variants.emplace_back((size_t)9, (size_t)2);

    ReadString pattern("AAAA");
    auto chromosomeLength = seqan::length(chromosome.data);
    auto startPos = 0;
    auto endPos = chromosomeLength;

    VariantIndex variantIndex(chromosome.variants, chromosomeLength - 1);

    SemiGlobalAligner aligner(variantIndex);

    std::vector<Alignment> alignments = aligner.align(chromosome, startPos, endPos, pattern, maxErrors);
    if (alignments.size() == 0) {
        std::cout << "No alignment with required error bound found!" << std::endl;
    } else {
        std::cout << "Best alignment is: " << alignments[0].getCigar() << std::endl;
        std::cout << "Alignment string is: " << alignments[0].getAlignmentString() << std::endl;
    }


//    ReferenceReader referenceReader(chromosomes, chromosomeFilter);
//    VariantsReader variantsReader(chromosomes, chromosomeFilter);
//    referenceReader.readReferences(referenceFiles);
//    variantsReader.readVariants(variantFiles);

//    chromosomes.push_back(chromosome);

    return 0;
}

