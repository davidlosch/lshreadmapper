#include <iostream>
#include <map>
#include "Stats.h"

void calculateStatistics(VariantStatisticsModule& vsm) {

    typedef size_t Pos;
    typedef size_t Count;
    typedef ptrdiff_t Length;
    typedef std::string Chromosome;

    std::vector<std::pair<Pos, Count> >& snps = vsm.snpCountPerPos;
    std::vector<std::pair<Pos, Count> >& indels = vsm.indelCountPerPos;
    std::vector<std::pair<Pos, Length> >& indelLength = vsm.indelLen;
    std::vector<std::pair<Chromosome, Pos> >& startPositions = vsm.chromosomeStartPositions;

    std::cout << "calculateStatistics" << std::endl;
    std::ofstream file("./stats.csv");

    std::vector<std::pair<Pos, Count>> snpAndIndels;
    snpAndIndels.reserve(snps.size() + indels.size());
    auto snpIt = snps.begin();
    auto snpEnd = snps.end();
    for (auto entry : indels) {
        while (snpIt != snpEnd && snpIt->first < entry.first) {
            snpAndIndels.push_back(*(snpIt++));
        }
        if (snpIt != snpEnd && snpIt->first == entry.first) {
            snpAndIndels.emplace_back(entry.first, entry.second + (snpIt++)->second);
        } else {
            snpAndIndels.push_back(entry);
        }
    }
    while (snpIt != snpEnd) {
        snpAndIndels.push_back(*(snpIt++));
    }
    snpAndIndels.shrink_to_fit();

    std::cout << "snps: " << snps.size() << "\n";
    std::cout << "indels: " << indels.size() << "\n";
    std::cout << "both: " << snpAndIndels.size() << "\n";
    std::cout << "last-pos: " << snpAndIndels.back().first << "\n";
    std::cout << std::endl;


    std::vector<int> q32_report(32);
    for (size_t i = 0; i < q32_report.size(); i++) {
        q32_report[i] = i+1;
    }
    std::vector<double> limits_report{std::numeric_limits<double>::max(), 65536, 16};
    writeMatrixAsCSV(
                histogram_numberOfCombinations(snps, q32_report, limits_report),
                file,
                "Combinations;SNPs only \nq=;limit=");

    std::vector<int> q16(16);
    for (size_t i = 0; i < q16.size(); i++) {
        q16[i] = i+1;
    }

    std::vector<int> q256(9);
    q256[0] = 1;
    for (size_t i = 1; i < q256.size(); i++) {
        q256[i] = q256[i-1] * 2;
    }
    std::vector<double> limits{16, 64, 256, 1024, 4096, 16384, 65536, std::numeric_limits<double>::max()};
    writeMatrixAsCSV(
                histogram_numberOfCombinations(snps, q16, limits),
                file,
                "Combinations;SNPs only \nq=;limit=");
    writeMatrixAsCSV(
                histogram_numberOfCombinations(indels, q16, limits),
                file,
                "Combinations;Indels only \nq=;limit=");
    writeMatrixAsCSV(
                histogram_numberOfCombinations(snps, q256, limits),
                file,
                "Combinations;SNPs only \nq=;limit=");
    writeMatrixAsCSV(
                histogram_numberOfCombinations(indels, q256, limits),
                file,
                "Combinations;Indels only \nq=;limit=");

    std::vector<size_t> distanceSteps;
    for (double d = 1; d < 10e7; d *= pow(10, 1./4.)) {
        distanceSteps.push_back(round(d));
    }
    writePairsAsCSV(
                convertMapToHistogram(
                    histogram_distanceBetweenVariants(snps),
                    distanceSteps),
                file,
                "Distance;SNPs only \nDistance;#");
    writePairsAsCSV(
                convertMapToHistogram(
                    histogram_distanceBetweenVariants(indels),
                    distanceSteps),
                file,
                "Distance;Indels only \nDistance;#");
    writePairsAsCSV(
                convertMapToHistogram(
                    histogram_distanceBetweenVariants(snpAndIndels),
                    distanceSteps),
                file,
                "Distance;SNPs and Indels \nDistance;#");

    writeHistogramAsCSV(
                histogram_numberOfVariantsInQgram(snps, 16),
                file,
                "q=16;SNPs only \nNumber of variants;#");
    writeHistogramAsCSV(
                histogram_numberOfVariantsInQgram(indels, 16),
                file,
                "q=16;Indels only \nNumber of variants;#");
    writeHistogramAsCSV(
                histogram_numberOfVariantsInQgram(snpAndIndels, 16),
                file,
                "q=16;SNPs and Indels \nNumber of variants;#");
    writeHistogramAsCSV(
                histogram_numberOfVariantsInQgram(snps, 201),
                file,
                "q=201;SNPs only \nNumber of variants;#");
    writeHistogramAsCSV(
                histogram_numberOfVariantsInQgram(indels, 201),
                file,
                "q=201;Indels only \nNumber of variants;#");
    writeHistogramAsCSV(
                histogram_numberOfVariantsInQgram(snpAndIndels, 201),
                file,
                "q=201;SNPs and Indels \nNumber of variants;#");

    writePairsAsCSV(
                insertAllXvaluesToMap(
                    histogram_indelLength(indelLength)),
                file,
                "IndelLength \nLength;#");

    writePairsAsCSV(
                insertAllXvaluesToMap(
                    histogram_numberOfVariantsPerPosition(snps)),
                file,
                "Number of variants per position;SNPs only \nNumber of variants;#");
    writePairsAsCSV(
                insertAllXvaluesToMap(
                    histogram_numberOfVariantsPerPosition(indels)),
                file,
                "Number of variants per position;Indels only \nNumber of variants;#");

    Pos pos = findOnePositionWithManyVariants(snpAndIndels, 201, 201);
    writePairsAsCSV(
                manyVariantsExampleWindow(snps, indels, pos, 201),
                file,
                "Many variants in one window (position: " + toString(pos) + ")\nSNPs;Indels");


    std::vector<Pos> allowedGaps = {0,1,2,3,4};
    std::vector<size_t> runLengthSteps;// = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000};
    runLengthSteps.push_back(1);
    for (double d = 1; d < 10000; d *= 1.1) {
        size_t value(d);
        if (value > runLengthSteps.back()) {
            runLengthSteps.push_back(d);
        }
    }

    writeMatrixAsCSV(
                histogram_noGapRunLength<double>(snps, allowedGaps, runLengthSteps),
                file,
                "Number of sequences with variants on all positions;SNPs only \nLength=;Allowed gap length=");
    writeMatrixAsCSV(
                histogram_noGapRunLength<double>(indels, allowedGaps, runLengthSteps),
                file,
                "Number of sequences with variants on all positions;Indels only \nLength=;Allowed gap length=");
    writeMatrixAsCSV(
                histogram_noGapRunLength<double>(snpAndIndels, allowedGaps, runLengthSteps),
                file,
                "Number of sequences with variants on all positions;SNPs and Indels \nLength=;Allowed gap length=");

    std::vector<std::pair<Chromosome, Pos> > chromosomeLength(startPositions.size() - 1);
    for (size_t i = 0; i < chromosomeLength.size(); i++) {
        chromosomeLength[i].first = startPositions[i].first;
        chromosomeLength[i].second = startPositions[i+1].second - startPositions[i].second - vsm.SPACE_BETWEEN_CHROMOSOMS;
    }
    writePairsAsCSV(
                chromosomeLength,
                file,
                "Chromosome length \nChromosome;Length");
//    writePairsAsCSV(
//                convertPosToChromosomePos(findPositionsWithManyVariants(snpAndIndels, 400, 200, 2000u), startPositions),
//                file,
//                "Positions with many variants (200 variants in 400 base pairs) \nChromosome;Position");
//    writePairsAsCSV(
//                convertPosToChromosomePos(findPositionsWithManyVariants(snpAndIndels, 100, 50, 2000u), startPositions),
//                file,
//                "Positions with many variants (50 variants in 100 base pairs) \nChromosome;Position");

    writeMatrixAsCSV(
                findPositionsWithManyVariants(snpAndIndels, startPositions, 200, 100),
                file,
                "Positions with many variants (100 variants in 200 base pairs) \nChromosome;Position;Length;Variants");
    writeMatrixAsCSV(
                findPositionsWithManyVariants(snpAndIndels, startPositions, 2000, 1000),
                file,
                "Positions with many variants (1000 variants in 2000 base pairs) \nChromosome;Position;Length;Variants")
    ; // so ein Pech aber auch, dass dieses Semikolon nicht mehr in die darüberliegende Zeile passte...

}


template <class T>
void showVector(const std::vector<T>& v) {
    using namespace std;
    for (size_t i = 0; i < v.size(); i++) {
        cout << i << ": " << v[i] << endl;
    }
}


void testStats() {
    std::map<size_t, size_t> m;
    m[3] = 1;
    m[5] = 1;
    m[6] = 2;
    m[7] = 5;
    m[10] = 1;

    typedef std::pair<size_t, size_t> Pair;
    std::vector<Pair> all;
    for (auto entry: m) {
        all.push_back(entry);
    }

    std::sort(all.begin(), all.end(), [](Pair a, Pair b){return a.first < b.first;});


    auto res1 = histogram_distanceBetweenVariants(all);
    writePairsAsCSV(res1, std::cout); // res[1..3] = { 2 1 1 }

    auto res2 = histogram_numberOfVariantsInQgram(all, 3);
    writeHistogramAsCSV(res2, std::cout); // res[0..3] = { 3 6 3 1 }

    auto res3 = histogram_numberOfCombinations<long long>(all, 5);
    writeHistogramAsCSV(res3, std::cout); // res[1..4] = { 21 44 ? ? }

    auto res4 = histogram_numberOfCombinations(all, 5, 4LL);
    writeHistogramAsCSV(res4, std::cout);

    std::vector<std::pair<unsigned, int> > indelLength;
    indelLength.emplace_back(3, 2);
    indelLength.emplace_back(5, 3);
    indelLength.emplace_back(6, 2);
    indelLength.emplace_back(7, 2);
    indelLength.emplace_back(9, 1);
    indelLength.emplace_back(15, -5);
    indelLength.emplace_back(16, 3);
    writePairsAsCSV(
                histogram_indelLength(indelLength),
                std::cout,
                "IndelLength \nLength;#"); // res[-5,1,2,3] = {1,1,3,2}
    writePairsAsCSV(
                insertAllXvaluesToMap(histogram_indelLength(indelLength)),
                std::cout,
                "IndelLength \nLength;#"); // all x values...

    writeMatrixAsCSV(
                histogram_noGapRunLength<size_t>(all, std::vector<size_t>{0,1,2}, std::vector<size_t>{1,2,3,4,5,6}),
                std::cout); // allowed=0: {2 0 1 0 0 0}, =1: {1 0 0 1 0 0}

    writeHistogramAsCSV(findPositionsWithManyVariants(all, 2, 2, (size_t)2), std::cout); // 5
    writeHistogramAsCSV(findPositionsWithManyVariants(all, 2, 2, (size_t)0), std::cout); // 5,6


    std::vector<std::pair<std::string, size_t> > startPositions;
    startPositions.emplace_back("A", 0);

    all.emplace_back(100,1);
    writeMatrixAsCSV(findPositionsWithManyVariants(all, startPositions, 2, 1), std::cout);

    all.emplace_back(101,1);
    writeMatrixAsCSV(findPositionsWithManyVariants(all, startPositions, 2, 2), std::cout);
}
