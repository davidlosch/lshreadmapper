#include "MapperModule.h"

MapperModule::MapperModule() {

}

void MapperModule::addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName) {
    using namespace seqan;

    addUsageLine(argParser, moduleName + " \\fIINDEX_FILE\\fP \\fIREADS.fastq\\fP [\\fI-p PAIRED_END_READS.fastq\\fP]");

    addSection(argParser, moduleName + " Options");
    ArgParseArgument indexFile(ArgParseArgument::INPUTFILE, "INDEX");
    addArgument(argParser, indexFile);
//    setValidValues(indexFile, "mvi");
    ArgParseArgument readsFile(ArgParseArgument::INPUTFILE, "READS FASTQ");
    addArgument(argParser, readsFile);
    setValidValues(readsFile, "fastq fq fastq.gz fq.gz");
    ArgParseOption pairedEndReadsFile("p", "paired-end", "Paired-end FASTQ read file.", ArgParseArgument::INPUTFILE);
    addOption(argParser, pairedEndReadsFile);
    setValidValues(pairedEndReadsFile, "fastq fq fastq.gz fq.gz");
}

void MapperModule::readOptions(seqan::ArgumentParser &argParser) {
    seqan::getArgumentValue(indexFilePath, argParser, 1);
    seqan::getArgumentValue(readsFilePath, argParser, 2);
    seqan::getOptionValue(pairedEndReadsFilePath, argParser, "p");
}

int MapperModule::run() {
    std::cout << "mapreads not yet implemented!" << std::endl;
    return -1;
}
