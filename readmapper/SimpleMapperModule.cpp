#include "SimpleMapperModule.h"
#include "DocumentSet.h"
#include "Timer.h"
#include "align/SemiGlobalAligner.h"

#include <cstdio>
#include <iostream>

SimpleMapperModule::SimpleMapperModule() {

}

void SimpleMapperModule::addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName) {
    using namespace seqan;

    addUsageLine(argParser, moduleName + " \\fIREFERENCE.fasta\\fP \\fIREADS.fastq\\fP");

//    addSection(argParser, moduleName + " Options");
    ArgParseArgument referenceFile(ArgParseArgument::INPUTFILE, "REFERENCE");
    addArgument(argParser, referenceFile);
    setValidValues(referenceFile, "fasta fa fasta.gz fa.gz");
    ArgParseArgument readsFile(ArgParseArgument::INPUTFILE, "READS FASTQ");
    addArgument(argParser, readsFile);
    setValidValues(readsFile, "fastq fq fastq.gz fq.gz");
}

void SimpleMapperModule::readOptions(seqan::ArgumentParser &argParser) {
    seqan::getArgumentValue(referenceFilePath, argParser, 1);
    seqan::getArgumentValue(readsFilePath, argParser, 2);
}

int SimpleMapperModule::run() {

    std::ofstream file_genome("./genom.fasta");
    std::mt19937 rand;
//    size_t size_genome = 1000*1000*30;
    size_t size_genome = 1000;
    seqan::Dna5String str;
    for (size_t i = 0; i < size_genome; i++) {
        str += rand() % 4;
    }
    seqan::writeRecord(file_genome, "0", str, seqan::Fasta());
    file_genome.flush();
    seqan::clear(str);


    seqan::SequenceStream refStream(referenceFilePath.c_str());



    seqan::CharString id;
    seqan::DnaString refGenome;
    seqan::readRecord(id, refGenome, refStream);

    DocumentSet documentSet;

    std::vector<Chromosome> chromosomes;

    Chromosome chromosome("testChromosome", refGenome, std::vector<Variant>());
//    chromosome.variants.emplace_back((size_t)15, (size_t)2);
    chromosomes.push_back(chromosome);

    auto chromosomeLength = seqan::length(chromosome.data);
    
    VariantIndex variantIndex(chromosome.variants, chromosomeLength - 1);

    SemiGlobalAligner aligner(variantIndex);

//    auto alignment = aligner.align(chromosome, startPos, endPos, pattern, maxErrors);
//    std::cout << "Best alignment is: " << alignment.getCigar() << std::endl;
    
    Timer timer;
    timer.startTimer();
    
    // //////////////////
    // Parameters
    
    const unsigned int maxErrors = 20;
    const unsigned int windowSize = 101;
    const unsigned int alignOffset = 0.5 * windowSize;
    documentSet.packageSize = 10000;
    
    // //////////////////

    documentSet.createHashFunctions();
    
    // Chromosome ID = 0
    documentSet.readFastaReferenceGenome(refGenome, windowSize, 0);

    std::ofstream file("./blub.fasta");
    size_t maximalErrorCountInTheReads = 4;
    for (int i = 0; i < 1000; i++) {
        seqan::Dna5String str;
        int off = rand() % (seqan::length(refGenome) - 100);
        for (int j = 0; j < 100; j++) {
            str += refGenome[j + off];
        }
        for (int j = 0; j < maximalErrorCountInTheReads; j++) {
            int r = rand() % 100;
            str[r] = unsigned(str[r]) ^ (rand() % 4);
        }

        seqan::writeRecord(file, std::to_string(i), str, seqan::Fasta());
    }
    file.flush();

    seqan::SequenceStream readsStream(readsFilePath.c_str());
    if (!isGood(refStream) || !isGood(readsStream)) {
        std::cerr << "ERROR: Could not open the file.\n";
        return 1;
    }
    
    auto matchedFunction = [&](std::unordered_set<MatchingPair> &matchingPairs){
        int i = 0;
        int max_i = matchingPairs.size();
        int alignedReads = 0;
        std::vector<int> alignedReadsWithError = std::vector<int>(maxErrors+1);
        std::vector<int> alignedReadsWithError2 = std::vector<int>(maxErrors+1, 0);
        std::unordered_map<int, int> errors;

        Timer alignTimer;
        alignTimer.startTimer();

        for (const MatchingPair &matchingPair : matchingPairs){
            unsigned int startPosition = documentSet.referenceDocumentPositions.at(matchingPair.referencePartIndex);
            startPosition = std::max(startPosition, alignOffset) - alignOffset;
            unsigned int endPosition = std::min((size_t)(startPosition + windowSize + 2*alignOffset), seqan::length(chromosome.data));
//            std::cout << "_:" << alignOffset << "," << documentSet.referenceDocumentPositions.at(matchingPair.referencePartIndex) << "," << startPosition << "," << endPosition
//                      << "," << windowSize << "," << alignOffset << "\n";
            auto &chromosome = chromosomes.at(matchingPair.referenceChromosome);
            std::vector<Alignment> alignments =
                aligner.align(chromosome,
                          startPosition,
                          endPosition,
                          documentSet.originalReadDocuments.at(matchingPair.readPartIndex),
                          maxErrors);

            if (alignments.size() != 0 && alignments[0].getErrorCount() > maximalErrorCountInTheReads) {
                std::cout << "Aligning pair " << ++i << " of " << max_i << ":" << std::endl;
                auto text = seqan::infix(chromosome.data, startPosition, endPosition);
                std::cout << "REF:  " << text << std::endl;
                std::cout << "READ: " << documentSet.originalReadDocuments.at(matchingPair.readPartIndex) << std::endl;
                for (Alignment a : alignments) {
                    std::cout << "ALN:  " << a.getAlignmentString() << std::endl;
                    std::cout << a.getCigar().c_str() << " (" << a.getErrorCount() << " errors)" << std::endl;
                }
                std::cout << "Read-Index: " << matchingPair.readPartIndex << std::endl;
                std::cout << std::endl;
            }
            if (alignments.size() > 0) {
                alignedReads++;
                alignedReadsWithError[alignments[0].getErrorCount()]++;
                if (errors.find(matchingPair.readPartIndex) == errors.end()) {
                    errors[matchingPair.readPartIndex] = alignments[0].getErrorCount();
                } else {
                    errors[matchingPair.readPartIndex] = std::min(errors[matchingPair.readPartIndex], alignments[0].getErrorCount());
                }

            }
        }
        int alignedReads2 = 0;
        //for (const std::pair<int,int>& err: errors) {
        for (unsigned int readIndex = 0; readIndex < documentSet.packageSize; readIndex++) {
            unsigned int errorCountOriginalRead = std::numeric_limits<unsigned int>::max();
            unsigned int errorCountReverseComplementRead = std::numeric_limits<unsigned int>::max();
            if (errors.find(readIndex) != errors.end()) {
                errorCountOriginalRead = errors[readIndex];
            }
            if (errors.find(readIndex + documentSet.packageSize) != errors.end()) {
                errorCountReverseComplementRead = errors[readIndex + documentSet.packageSize];
            }
            unsigned int errorCount = std::min(errorCountOriginalRead, errorCountReverseComplementRead);
            if (errorCount != std::numeric_limits<unsigned>::max()) {
                alignedReadsWithError2[errorCount]++;
                alignedReads2++;

                if (errorCount > maximalErrorCountInTheReads) {
                    std::cout << "Too many errors: Read-Index = " << readIndex << std::endl;
                }
            }
        }

        timer.stopTimer();
        std::cout << "Aligned " << alignedReads2 << " of " << max_i << " reads." << std::endl;
        for (int j = 0; j <= maxErrors; j++) {
            std::cout << "Alignment with " << j << " errors: " << alignedReadsWithError2[j] << std::endl;
        }
        timer.printDuration("Alignment time");
        std::cout << std::endl;
    };
    auto unmatchedFunction = [&](std::vector<unsigned int> &unmatchedReads){
        /*
        std::cout << "Unmatched Reads: " << unmatchedReads.size() << std::endl;
        for (unsigned int readId : unmatchedReads) {
            std::cout << readId << std::endl;
            std::cout << documentSet.originalReadDocuments.at(readId) << std::endl;
        }
        */
    };

    documentSet.readFastQReads(readsStream, timer, matchedFunction, unmatchedFunction);

    timer.printAndRestartTimer("Q-Gram and signature creation");


//    documentSet.findSimilarDocuments();

    timer.printAndRestartTimer("Similarity search");

    // TODO: write output!
#if 0
    // 4) Write out Sam file.
    std::ofstream samFile(argv[3], std::ios_base::out);
    write(samFile, fragStore, Sam());
#endif

    return 0;
}
