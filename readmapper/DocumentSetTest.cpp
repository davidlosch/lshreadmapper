#include <iostream>
#include "DocumentSet.h"
#include "DocumentSetTest.h"
#include <thread>

seqan::Iupac toIupac(int t, int a, int c, int g) {
    return seqan::Iupac(t + (a << 1) + (c << 2) + (g << 3));
}


bool test_findCombinationsOfQgram() {

    std::vector<seqan::Iupac> iupacString;
    iupacString.push_back(toIupac(0, 0, 1, 1));
    iupacString.push_back(toIupac(0, 0, 0, 1));
    iupacString.push_back(toIupac(0, 0, 1, 0));
    iupacString.push_back(toIupac(0, 1, 0, 0));
    iupacString.push_back(toIupac(1, 1, 1, 1));
    iupacString.push_back(toIupac(1, 0, 0, 0));
    iupacString.push_back(toIupac(0, 1, 1, 1));
    iupacString.push_back(toIupac(0, 0, 1, 1));

    int count = 0;

    typedef typename std::vector<seqan::Dna>::iterator OpIt;
    auto func = [&count](OpIt beginIt, OpIt endIt) {
        count++;
//        for (; beginIt != endIt; ++beginIt) {
//            std::cout << seqan::Dna(*beginIt);
//        }
//        std::cout << std::endl;
    };


    std::vector<seqan::Dna> space1(iupacString.size());
    std::vector<char> space2(iupacString.size() + 1);
    findCombinationsOfQGram(iupacString, 0, iupacString.size(), space1, space2, func);

    return count == 48;
}

