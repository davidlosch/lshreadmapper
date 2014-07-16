#ifndef PG583_MAPPER_MODULE_H
#define PG583_MAPPER_MODULE_H

#include "ApplicationModule.h"

#include <string>

class MapperModule : public ApplicationModule {
private:
    std::string indexFilePath;
    std::string readsFilePath;
    std::string pairedEndReadsFilePath;
public:
    MapperModule();
    void addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName);
    void readOptions(seqan::ArgumentParser &argParser);
    int run();
};

#endif // PG583_MAPPER_MODULE_H
