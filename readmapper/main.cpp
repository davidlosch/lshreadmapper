#include "module/AlignmentTestModule.h"
#include "module/VariantStatisticsModule.h"
#include "module/MapperModule.h"
#include "module/IndexerModule.h"

#include <seqan/arg_parse.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

int main(int argc, const char *argv[]) {
    IndexerModule indexer;
    MapperModule mapper;
    VariantStatisticsModule vcfStats;
    AlignmentTestModule alignTest;

    std::vector<std::reference_wrapper<ApplicationModule>> modules{indexer, mapper, vcfStats, alignTest};
    std::vector<std::string> moduleNames;
    for (ApplicationModule &module : modules) {
        moduleNames.push_back(module.getName());
    }

    std::string applicationName = "readmapper";
    seqan::ArgumentParser argParser(applicationName);
    seqan::setShortDescription(argParser, "A variant tolerant read mapper and aligner.");
    seqan::addArgument(argParser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING));
    auto &moduleNameArg = seqan::getArgument(argParser, 0);
    seqan::setValidValues(moduleNameArg, moduleNames);

    for (ApplicationModule &module : modules) {
        seqan::addUsageLine(argParser, module.getName() + " [--help | \\fImodule-options\\fP]");
    }
//    seqan::addDescription(argParser, "description");
    std::ostringstream out, err;
    seqan::ArgumentParser::ParseResult parseRes = seqan::parse(argParser, argc, argv, out, err);
    if (!seqan::hasValue(moduleNameArg)) {
        std::cout << out.str();
        std::cerr << err.str();
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    ApplicationModule &module = modules[std::find(moduleNames.begin(), moduleNames.end(),
                                                   seqan::getArgumentValue(moduleNameArg)) - moduleNames.begin()];
    seqan::ArgumentParser moduleArgParser(applicationName);
    seqan::addArgument(moduleArgParser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, module.getName()));
    module.addOptions(moduleArgParser);
    seqan::addUsageLine(moduleArgParser, module.generateUsageLine(moduleArgParser));

    parseRes = seqan::parse(moduleArgParser, argc, argv);

    if (parseRes != seqan::ArgumentParser::PARSE_OK) {
        // outputs predefined functions like 'help', 'version' etc. if not PARSE_ERROR
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    module.readOptions(moduleArgParser);

    return module.run();
}
