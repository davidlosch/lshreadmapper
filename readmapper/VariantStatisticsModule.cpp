#include "VariantStatisticsModule.h"
#include "Stats.h"

#include <seqan/vcf_io.h>
#include <seqan/basic.h>

#include <functional>
#include <iostream>
#include <vector>
#include <utility>
#include <cstdio>
#include <cstdint>

VariantStatisticsModule::VariantStatisticsModule() {

}

void VariantStatisticsModule::addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName) {
    using namespace seqan;

    addUsageLine(argParser, moduleName + " [\\fIOPTIONS\\fP] \\fISORTED_VARIANTS.vcf\\fP");

    addSection(argParser, moduleName + " Options");
    ArgParseArgument variantsFile(ArgParseArgument::INPUTFILE, "SORTED VARIANTS VCF");
    addArgument(argParser, variantsFile);
    setValidValues(variantsFile, "vcf vcf.gz");
}

void VariantStatisticsModule::readOptions(seqan::ArgumentParser &argParser) {
    seqan::getArgumentValue(variantsFilePath, argParser, 1);
}

int VariantStatisticsModule::run() {
    snpCountPerPos.reserve(65000000);
    indelCountPerPos.reserve(10000000);
    indelLen.reserve(10000000);

    if(readVariantFile() != 0) {
        return 1;
    }
    snpCountPerPos.shrink_to_fit();
    indelCountPerPos.shrink_to_fit();
    indelLen.shrink_to_fit();

    createStatistics();
    return 0;
}

int VariantStatisticsModule::readVariantFile() {
    std::hash<std::string> hash_fn;

    size_t totalPosCount = 0;
    size_t totalSNPCount = 0;
    size_t totalInsertionCount = 0;
    size_t totalDeletionCount = 0;

    size_t prevPos = 0;
    size_t chromosomOffset = 0;
    size_t prevChromosomHash = hash_fn("");
    size_t curSNPCount = 0;
    size_t curIndelCount = 0;

    std::ifstream file(variantsFilePath.c_str());
    if (!file) {
        std::cerr << "ERROR: Could not open " << variantsFilePath << std::endl;
        return 1;
    }
    chromosomeStartPositions.emplace_back("1", 0); // first chromosome starts at position 0
    while (file.good()) {
        std::string lineStr;
        std::getline(file, lineStr);
        if (lineStr.empty() || lineStr.at(0) == '#') {
            continue;
        }
        std::stringstream line(lineStr);
        std::string chromosom;
        std::getline(line, chromosom, '\t');
        std::string posStr;
        std::getline(line, posStr, '\t');
        std::string id;
        std::getline(line, id, '\t');
        std::string ref;
        std::getline(line, ref, '\t');
        std::string altStr;
        std::getline(line, altStr, '\t');

        if(chromosom.empty() || posStr.empty() || ref.empty() || altStr.empty()) {
            continue;
        }

        std::size_t chromosomHash = hash_fn(chromosom);
        if (totalPosCount != 0 && (chromosomHash != prevChromosomHash || totalPosCount % 100000 == 0)) {
            std::cout << "positions: " << std::setw(8) << totalPosCount <<
                         "  variants: " << std::setw(8) << (totalSNPCount + totalInsertionCount + totalDeletionCount) <<
                         "  SNPs: " << std::setw(8) << totalSNPCount <<
                         "  Indels: " << std::setw(8) << (totalInsertionCount + totalDeletionCount) <<
                         "  Insertions: " << std::setw(8) << totalInsertionCount <<
                         "  Deletions: " << std::setw(8) << totalDeletionCount <<
                         "  chromosome: " << std::setw(2) << chromosom <<
                         "  pos-in-chromosome: " << std::setw(9) << posStr <<
                         "  prevPos: " << std::setw(9) << prevPos <<
                         std::endl;
        }

        if (chromosomHash != prevChromosomHash) {
            // space chromosomes apart
            if (totalPosCount != 0) {
                chromosomOffset = prevPos + SPACE_BETWEEN_CHROMOSOMS;
                chromosomeStartPositions.emplace_back(chromosom, chromosomOffset);
            }
            prevChromosomHash = chromosomHash;
        }
        size_t pos = chromosomOffset + (size_t) std::stoul(posStr);
        if (pos == prevPos) {
            if(curSNPCount != 0) {
                snpCountPerPos.pop_back();
            }
            if(curIndelCount != 0) {
                indelCountPerPos.pop_back();
            }
        } else {
            curSNPCount = 0;
            curIndelCount = 0;
            ++totalPosCount;
        }

        size_t refSize = ref.size();
        std::stringstream alts(altStr);
        std::string alt;
        while (std::getline(alts, alt, ',')) {
            if (!alt.empty()) {
                size_t altSize = alt.size();
                if (altSize != refSize) {
                    ++curIndelCount;
                    if(altSize > refSize) {
                        ++totalInsertionCount;
                        indelLen.emplace_back(pos, altSize - refSize);
//                        indelLen.emplace_back(pos, altSize);
                    } else {
                        ++totalDeletionCount;
                        indelLen.emplace_back(pos, altSize - refSize);
//                        indelLen.emplace_back(pos, -altSize);
                    }
                } else {
                    ++curSNPCount;
                    ++totalSNPCount;
                }
            }
        }
        if (curSNPCount != 0) {
            if (snpCountPerPos.size() == 0 || pos > snpCountPerPos.back().first) {
                snpCountPerPos.emplace_back(pos, curSNPCount);
            } else if (pos == snpCountPerPos.back().first) {
                snpCountPerPos.back().second += curSNPCount;
            } else {
                std::cout << "input is not sorted: " << pos << "," << snpCountPerPos.back().first << "\n";
            }
        }
        if (curIndelCount != 0) {
            if (indelCountPerPos.size() == 0 || pos > indelCountPerPos.back().first) {
                indelCountPerPos.emplace_back(pos, curIndelCount);
            } else if (pos == indelCountPerPos.back().first) {
                indelCountPerPos.back().second += curIndelCount;
            } else {
                std::cout << "input is not sorted: " << pos << "," << indelCountPerPos.back().first << "\n";
            }
        }
        prevPos = pos;
        //if (pos > 1000000000) break;
    }
    std::cout << "----------------" << std::endl;
    std::cout << "positions: " << std::setw(8) << totalPosCount <<
                 "  variants: " << std::setw(8) << (totalSNPCount + totalInsertionCount + totalDeletionCount) <<
                 "  SNPs: " << std::setw(8) << totalSNPCount <<
                 "  Indels: " << std::setw(8) << (totalInsertionCount + totalDeletionCount) <<
                 "  Insertions: " << std::setw(8) << totalInsertionCount <<
                 "  Deletions: " << std::setw(8) << totalDeletionCount <<
                 "  last-pos: " << std::setw(9) << prevPos <<
                 std::endl;
    return 0;
//    // <seqan/vcf_io> does not work currently, do not know why...
//    using namespace seqan;
//    VcfStream vcfIn(variantsFilePath.c_str());
//    if (!isGood(vcfIn)) {
//        std::cerr << "ERROR: Could not open " << variantsFilePath << std::endl;
//        return 1;
//    }
//
//    VcfRecord record;
//    while (!atEnd(vcfIn)) {
//        if (readRecord(record, vcfIn) != 0) {
//            std::cerr << "ERROR: Problem reading from " << variantsFilePath << std::endl;
//            return 1;
//        }
//    }
}

void VariantStatisticsModule::createStatistics() {
    calculateStatistics(*this);
}
