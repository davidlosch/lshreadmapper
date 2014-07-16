#ifndef LSH_H
#define LSH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "Chromosome.h"



class LSH {
public:
    /**
     * @brief Represents a window of the genome
     * @details This struct stores a pointer to the respective chromosome and the begin position of the window (aka document)
     */
    struct Reference {
        Chromosome* chromosome;
        unsigned beginPosition;
    };

    typedef unsigned ReferenceID;
    typedef unsigned BandHashValue;

    size_t windowSize;
    size_t windowOverlappingOffset;
    size_t qGramSize;
    size_t signatureLength; // In other words: hash function count
    size_t bandSize; // For Locality-Sensitive Hashing
    size_t limit; // The number of combinations that is at most added
    size_t limit_skip; // If the number of combinations in a q-gram is bigger than this value, no combination is added

    std::minstd_rand rngLimit; // this random generator is used to select some q-grams if there are more qgrams than "limit". The generator should be fast, so don't use std::mt19937

    std::vector<unsigned> hashFunctions;
    std::vector<Reference> references; //
    std::unordered_map<BandHashValue, std::vector<ReferenceID>> bandHashTable;

    LSH();
    LSH(unsigned long long seed);



    /**
     * @brief Saves the LSH-index into a output-stream.
     * @details All attributes are saved.
     */
    void saveIndex(std::ostream&);

    void loadIndex(std::istream&);


    template <class Random>
    void createHashFunctions(Random& random);

    /**
     * @brief Adds all chromosomes between begin_it and end_it.
     */
    template <class Iterator>
    void addGenome(Iterator begin_it, Iterator end_it);


    /**
     * @brief
     * @details Creates overlapping windows of the chromosome (using the attributes windowSize and windowOverlappingOffset).
     *      Each window is stored as a Reference object and is added to the hash table via addReference().
     * @param chromosome
     */
    void addChromosome(Chromosome& chromosome);

    /**
     * @brief Adds one window of the reference genome.
     * @details This function adds the reference to the attribute references, calculates the q-grams via createQGrams_iupac,
     *      calculates their signatures and adds the corresponding band values to the bandHashTable
     * @param reference a window of the reference genome
     * @return the reference ID of the added reference.
     */
    ReferenceID addReference(Chromosome& chromosome, size_t beginPosition);

    /**
     * @brief creates the q-grams of a DnaString or Dna5String
     * @param window_begin
     * @param window_end
     * @param out_qGramHashValues After execution of this function, this variable contains the q-grams (respectively their hash values)
     */
    template <class Iterator>
    void createQGrams(Iterator window_begin, Iterator window_end, std::vector<unsigned>& out_qGramHashValues );

    /**
     * @brief creates the q-grams of a IupacString
     * @param window_begin
     * @param window_end
     * @param out_qGramHashValues After execution of this function, this variable contains the q-grams (respectively their hash values)
     */
    template <class Iterator>
    void createQGrams_iupac(Iterator window_begin, Iterator window_end, std::unordered_set<unsigned>& out_qGramHashValues);

    template <class Iterator>
    void createSignatures(Iterator qgrams_begin, Iterator qgrams_end, std::vector<unsigned>& out_signatures);

    template <class Iterator>
    void createBandValues(Iterator signature_begin, Iterator signature_end, std::vector<BandHashValue>& out_bandValues);

    template <class Iterator>
    void addBandValues(Iterator signature_begin, Iterator signature_end, ReferenceID& referenceID);

    /**
     * @brief Finds the read and its reverse complement in the genome.
     * @param out_referenceIDs Positions where the read was found.
     * @param out_referenceIDs_revComplement Positions where the reverse complement of the read was found.     */
    void findRead(ReadString read,
                  std::vector<ReferenceID>& out_referenceIDs,
                  std::vector<ReferenceID>& out_referenceIDs_revComplement);

    /**
     * @brief Finds the given string in the genome. The reverse complement is not taken into account
     * @param read_begin
     * @param read_end
     * @param out_referenceIDs Positions where the read was found.
     * @param out_referenceIDs_revComplement Positions where the reverse complement of the read was found.     */
    template <class Iterator>
    void findReadDirected(Iterator read_begin, Iterator read_end, std::vector<ReferenceID>& out_referenceIDs);


    template <class Iterator> // an iterator pointing to seqan::Dna or seqan::Dna5
    unsigned qGramHashFucntion(Iterator qgram_begin, Iterator qgram_end);

    template <class Iterator>
    BandHashValue bandHashFunction(Iterator signature_begin, Iterator signature_end);

    /** Find all variant combinations of the substring "iupacString[beginPos...beginPos+qGramSize]"
     *  and execute the function operation on the founds.
     *  @param curQGram arbitrary space of length qGramSize
     *  @param curVariant arbitrary space of length qGramSize+1
     *  @param operation a object that overloads operator()(Iterator, Iterator) with Iterator = std::vector<seqan::DNA>::iterator
     */
    template <class Iterator, class Function>
    void findCombinationsOfQGram(Iterator iupacQgram_begin, Function& operation);
private:
    std::vector<unsigned> tmp; // this vector is used as output parameter to avoid repeated memory allocations.

    //
    std::vector<seqan::Dna> curQGram;
    std::vector<char> curVariant;
};

// =====================================================================================
// ===== I M P L E M E N T A T I O N S
// =====================================================================================

template <class Iterator>
void LSH::addGenome(Iterator begin_it, Iterator end_it) {
    for (; begin_it != end_it; ++begin_it) {
        addChromosome(*begin_it);
    }
}

template <class Iterator>
void LSH::createQGrams(Iterator window_begin, Iterator window_end, std::vector<unsigned>& out_qGramHashes) {
    out_qGramHashes.clear();
    out_qGramHashes.reserve(window_end - window_begin - qGramSize + 1);
    for (auto qgram_end = window_begin + qGramSize; qgram_end <= window_end; ++window_begin, ++qgram_end) {
        // window_begin is the start position of the current q-gram.
        out_qGramHashes.push_back(qGramHashFucntion(window_begin, qgram_end));
    }
    // Remove duplicates.
    std::sort(out_qGramHashes.begin(), out_qGramHashes.end());
    out_qGramHashes.erase(std::unique(out_qGramHashes.begin(), out_qGramHashes.end()), out_qGramHashes.end());
}

template <class Iterator>
void LSH::createQGrams_iupac(Iterator window_begin, Iterator window_end, std::unordered_set<unsigned>& out_qGramHashValues) {
    // number of variants for a iupac symbol
    const static auto nv = [](unsigned x) {return (x & 1)
                                                + ((x & 2) >> 1)
                                                + ((x & 4) >> 2)
                                                + ((x & 8) >> 3);};
    const static unsigned long long NUMBER_OF_VARIANTS[16] = {nv(0), nv(1), nv(2), nv(3), nv(4), nv(5), nv(6), nv(7),
                                                              nv(8), nv(9), nv(10),nv(11),nv(12),nv(13),nv(14),nv(15)};

    curQGram.resize(qGramSize);
    curVariant.resize(qGramSize+1);

    out_qGramHashValues.clear();
    out_qGramHashValues.reserve(3 * qGramSize);
    typedef typename std::vector<seqan::Dna>::iterator OpIt;

    auto qgram_end = window_begin + qGramSize;
    unsigned long long combinations = 1; // the current number of combinations
    for (auto it = window_begin; it < qgram_end; ++it) {
        combinations *= NUMBER_OF_VARIANTS[unsigned(seqan::Iupac(*it))];
    }

    // the normal operation to add all combinations
    auto operation = [&](OpIt beginIt, OpIt endIt) {
        out_qGramHashValues.insert(qGramHashFucntion(beginIt, endIt));
    };

    // the operation that adds on average only "limit" combinations
    const double rngMax = rngLimit.max();
    const double limitAsFloat = limit;
    auto operation_limit = [&out_qGramHashValues, this, rngMax, &combinations, limitAsFloat](OpIt beginIt, OpIt endIt) {
        double r = this->rngLimit();
        if (r / rngMax < limitAsFloat / combinations) {
            out_qGramHashValues.insert(qGramHashFucntion(beginIt, endIt));
        }
    };

    auto lastQgram_begin = window_end - qGramSize;
    for (auto curQgram_begin = window_begin; curQgram_begin <= lastQgram_begin; ++curQgram_begin) {
        if (combinations <= limit) {
            // add all combinations
            findCombinationsOfQGram(curQgram_begin, operation);
        } else if (combinations <= limit_skip) {
            // add a random set of (on average) "limit" combinations
            findCombinationsOfQGram(curQgram_begin, operation_limit);
        } else {
            // skip this q-gram (because there are too many combinations)
        }

        // update the number of combinations for the next q-gram
        if (curQgram_begin < lastQgram_begin) { // true, if this is not the last iteration
            combinations /= NUMBER_OF_VARIANTS[unsigned(seqan::Iupac(*curQgram_begin))]; // there is no remainder for this division!
            combinations *= NUMBER_OF_VARIANTS[unsigned(seqan::Iupac(*(curQgram_begin + qGramSize)))];
        }
    }
}

template <class Iterator>
void LSH::createSignatures(Iterator qgrams_begin, Iterator qgrams_end, std::vector<unsigned>& out_signatures) {
    out_signatures.resize(signatureLength, std::numeric_limits<unsigned int>::max());
    for (; qgrams_begin != qgrams_end; ++qgrams_begin) {
        unsigned qGramHash = *qgrams_begin; //hashFunction(qgram);
        for (size_t hash = 0; hash < signatureLength; hash++) {
            // Create a normal hash value and XOR it with one of the randomly
            // generated integers. This will serve well enough as
            // "another" hash function
            unsigned int hashValue = qGramHash ^ hashFunctions[hash];
            out_signatures[hash] = hashValue < out_signatures[hash] ? hashValue : out_signatures[hash];
        }
    }
}

template <class Iterator>
void LSH::createBandValues(Iterator sig_begin, Iterator sig_end, std::vector<BandHashValue>& out_bandValues) {
    out_bandValues.clear();
    out_bandValues.reserve((sig_end - sig_begin) / 3);
    for (; sig_begin < sig_end; sig_begin = sig_begin + bandSize) {
        out_bandValues.push_back(bandHashFunction(sig_begin, sig_begin + bandSize));
    }
}

template <class Iterator>
void LSH::addBandValues(Iterator sig_begin, Iterator sig_end, LSH::ReferenceID& referenceID) {
    auto band_end = sig_begin + bandSize;
    for (; band_end < sig_end;) {
        bandHashTable[bandHashFunction(sig_begin, band_end)].push_back(referenceID);
        sig_begin = band_end;
        band_end = sig_begin + bandSize;
    }

}

template <class Iterator>
void LSH::findReadDirected(Iterator read_begin, Iterator read_end, std::vector<LSH::ReferenceID>& out_referenceIDs) {
    std::vector<unsigned> qgrams;
    createQGrams(read_begin, read_end, qgrams);
    std::vector<unsigned> signatures;
    createSignatures(qgrams.begin(), qgrams.end(), signatures);
    createBandValues(signatures.begin(), signatures.end(), out_referenceIDs);
}


template <class Iterator> // an iterator pointing to seqan::Dna or seqan::Dna5
unsigned LSH::qGramHashFucntion(Iterator begin, Iterator end) {
    static_assert(std::is_same<typename std::remove_reference<decltype(*begin)>::type, seqan::Dna>::value ||
                  std::is_same<typename std::remove_reference<decltype(*begin)>::type, seqan::Dna5>::value,
                  "Iterator is pointing to the wrong type");
    unsigned result = *begin++;
    for (; begin != end; ++begin) {
        result = result * 5 + unsigned(*begin);
    }
    return result;
}

template <class Iterator>
LSH::BandHashValue LSH::bandHashFunction(Iterator begin, Iterator end) {
    unsigned result = *begin++;
    for (; begin != end; ++begin) {
        result = result * 33767 + unsigned(*begin); // 33767 is a prime number
    }
    return result;
}


template <class Iterator, class Function>
void LSH::findCombinationsOfQGram(Iterator iupacQgram_begin, Function& operation) {
    int curIndex = 0; // current index in the current q-gram
    curVariant[0] = 0; // begin with variant id 0 at position 0

    const int qGramSizeMinusOne = (int) qGramSize - 1;

    while (true) {
        while (curVariant[curIndex] == 4) {
            // Alle Varianten wurden für das aktuelle Zeichen abgearbeitet
            curIndex--;
            if (curIndex < 0) {
                // Alle Varianten des aktuellen Q-Gramms wurden berücksichtigt -> break und nächstes Q-Gramm
                return;
            }
        }

        // extract the current iupac-base
        seqan::Iupac iupacSymbol = *(iupacQgram_begin + curIndex);
        seqan::Iupac dnaChar = iupacSymbol.value & (1 << curVariant[curIndex]);
        curVariant[curIndex]++; // next variant for the current position
        if (dnaChar.value != 0) {
            // an der Position "curIndex" im Q-Gramm kann die Base mit ID "curVariant[curIndex]" stehen
            curQGram[curIndex] = dnaChar;
            if (curIndex == qGramSizeMinusOne) {
                // we are at the end of the q-gram, so add the current q-gram to the document-q-gram-set
                operation(curQGram.begin(), curQGram.end());
           } else {
                curIndex++; // next position
                curVariant[curIndex] = 0; // current variant at the new position is set to zero
            }
        }
    }

}



#endif // LSH_H
