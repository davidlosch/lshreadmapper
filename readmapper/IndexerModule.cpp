#include "IndexerModule.h"

IndexerModule::IndexerModule() {

}

void IndexerModule::addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName) {
    using namespace seqan;

    addUsageLine(argParser, moduleName + " [\\fIOPTIONS\\fP] \\fIREFERENCE.fasta\\fP \\fIVARIANTS.vcf\\fP");

    addSection(argParser, moduleName + " Options");
    ArgParseArgument referenceFile(ArgParseArgument::INPUTFILE, "REFERNCE FASTA");
    addArgument(argParser, referenceFile);
    setValidValues(referenceFile, "fasta fa fasta.gz fa.gz");
    ArgParseArgument variantsFile(ArgParseArgument::INPUTFILE, "VARIANTS VCF");
    addArgument(argParser, variantsFile);
    setValidValues(variantsFile, "vcf vcf.gz"); // "vcf vcf.gz bcf" ???
    addOption(argParser, ArgParseOption("i", "index", "Write index to FILE.", ArgParseArgument::OUTPUTFILE));

    ArgParseOption qGramLength("q", "q-gram-length", "Use q-grams of length N.", ArgParseArgument::INTEGER, "N");
    addOption(argParser, qGramLength);
    setMinValue(qGramLength, "0");
    setMaxValue(qGramLength, "50");

    ArgParseOption editDistance("e", "edit-distance", "Permit an edit distance of up to N.",
                                ArgParseArgument::INTEGER, "N");
    addOption(argParser, editDistance);
    setMinValue(editDistance, "0");
    setMaxValue(editDistance, "10");
}

void IndexerModule::readOptions(seqan::ArgumentParser &argParser) {
    seqan::getArgumentValue(referenceFilePath, argParser, 1);
    seqan::getArgumentValue(variantsFilePath, argParser, 2);
    seqan::getOptionValue(qGramLength, argParser, "q");
    seqan::getOptionValue(editDistance, argParser, "e");
}

int IndexerModule::run() {
    std::cout << "creatindex not yet implemented!" << std::endl;
    return -1;
}
