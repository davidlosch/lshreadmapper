#ifndef PG583_VARIANT_STATISTICS_MODULE_H
#define PG583_VARIANT_STATISTICS_MODULE_H

#include "ApplicationModule.h"

#include <cstddef>

class VariantStatisticsModule : public ApplicationModule {
private:
    std::string variantsFilePath;

    int readVariantFile();
    void createStatistics();
public:
    std::vector<std::pair<size_t, size_t>> snpCountPerPos;
    std::vector<std::pair<size_t, size_t>> indelCountPerPos;
    std::vector<std::pair<size_t, ptrdiff_t>> indelLen;
    std::vector<std::pair<std::string, size_t>> chromosomeStartPositions;

    const static size_t SPACE_BETWEEN_CHROMOSOMS = 10000;

    VariantStatisticsModule();
    void addOptions(seqan::ArgumentParser &argParser, const std::string &moduleName);
    void readOptions(seqan::ArgumentParser &argParser);
    int run();
};

#endif // PG583_VARIANT_STATISTICS_MODULE_H
