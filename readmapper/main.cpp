#include "AlignmentTestModule.h"
#include "VariantStatisticsModule.h"
#include "SimpleMapperModule.h"
#include "MapperModule.h"
#include "IndexerModule.h"
#include "Stats.h"
#include "DocumentSetTest.h"

#include <seqan/arg_parse.h>

#include <sstream>
#include <string>
#include <map>

// mapreads-simple ../genome/phix_174/referenz/phix.fasta ../genome/phix_174/reads/phix.10k.fastq.gz

//#include "Test.h"
int main(int argc, const char *argv[]) {
//    Test::testJQ(); return 0;
    std::cout << "test_findCombinationsOfQgram(): " << test_findCombinationsOfQgram() << std::endl;
    //testStats(); return 0;
    seqan::ArgumentParser argParser("readmapper");

    IndexerModule indexer;
    MapperModule mapper;
    SimpleMapperModule simpleMapper;
    VariantStatisticsModule vcfStat;
    AlignmentTestModule alignTest;

    std::map<std::string, ApplicationModule*> modules;
    modules["createindex"] = &indexer;
    modules["mapreads"] = &mapper;
    modules["mapreads-simple"] = &simpleMapper;
    modules["vcfstat"] = &vcfStat;
    modules["testalign"] = &alignTest;

    seqan::ArgParseArgument mode(seqan::ArgParseArgument::STRING, "mode");
    seqan::addArgument(argParser, mode);

    std::stringstream ss;
    for(auto module: modules) {
        ss << module.first << " ";
    }
    seqan::setValidValues(mode, ss.str());

    if (argc < 2 || modules.count(argv[1]) == 0) {
        for(auto module: modules) {
            module.second->addOptions(argParser, module.first);
        }
        return seqan::parse(argParser, argc, argv);
    }
    const std::string &moduleName = modules.find(argv[1])->first;
    ApplicationModule &module = *modules.find(argv[1])->second;
    module.addOptions(argParser, moduleName);

    seqan::ArgumentParser::ParseResult argParseRes = seqan::parse(argParser, argc, argv);

    if (argParseRes != seqan::ArgumentParser::PARSE_OK) {
        // outputs predefined functions like 'help', 'version' etc. if not PARSE_ERROR
        return argParseRes != seqan::ArgumentParser::PARSE_ERROR;
    }

    module.readOptions(argParser);

    return module.run();
}
