#ifndef PG583_SIMPLE_MAPPER_MODULE_H
#define PG583_SIMPLE_MAPPER_MODULE_H

#include "ApplicationModule.h"

#include <string>

class SimpleMapperModule : public ApplicationModule {
private:
    std::string referenceFilePath;
    std::string readsFilePath;
public:
    SimpleMapperModule();
    void addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName);
    void readOptions(seqan::ArgumentParser &argParser);
    int run();
};

#endif // PG583_SIMPLE_MAPPER_MODULE_H
