#ifndef PG583_APPLICATION_MODULE_H
#define PG583_APPLICATION_MODULE_H

#include <seqan/arg_parse.h>

#include <string>

class ApplicationModule {
public:
    virtual void addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName) = 0;
    virtual void readOptions(seqan::ArgumentParser &argParser) = 0;
    virtual int run() = 0;
};

#endif // PG583_APPLICATION_MODULE_H
