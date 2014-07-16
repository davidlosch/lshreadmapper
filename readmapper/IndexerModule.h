#ifndef PG583_INDEXER_MODULE_H
#define PG583_INDEXER_MODULE_H

#include "ApplicationModule.h"

#include <string>

class IndexerModule : public ApplicationModule {
private:
    std::string referenceFilePath;
    std::string variantsFilePath;
    int qGramLength;
    int editDistance;
public:
    IndexerModule();
    void addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName);
    void readOptions(seqan::ArgumentParser &argParser);
    int run();
};

#endif // PG583_INDEXER_MODULE_H
