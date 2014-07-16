#ifndef PG583_ALIGNMENTTESTMODULE_H
#define PG583_ALIGNMENTTESTMODULE_H

#include "ApplicationModule.h"
#include "align/Variant.h"
#include "types.h"
#include "Logger.h"

#include <seqan/vcf_io.h>

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

class AlignmentTestModule : public ApplicationModule {
    std::unordered_set<std::string> chromosomeFilter;
    std::vector<std::string> variantFiles;
    std::vector<std::string> referenceFiles;
    int maxErrors;
public:
    AlignmentTestModule();
    void addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName);
    void readOptions(seqan::ArgumentParser &argParser);
    int run();
};

#endif // PG583_ALIGNMENTTESTMODULE_H
