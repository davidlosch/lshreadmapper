#ifndef PG583_DOCUMENTSET_H
#define PG583_DOCUMENTSET_H

#include "Timer.h"
#include "IupacString.h"

#include "seqan/sequence.h"
#include "seqan/seq_io.h"

#include <vector>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <random>

enum QGramType {
    Q_GRAM_WORD, Q_GRAM_CHAR
};

enum HashBucketEntryType {
    REFERENCE_PART, READ
};

struct HashBucketEntry {
    unsigned int documentID;
    unsigned int chromosome;
//    unsigned int bandID;
    bool operator==(const HashBucketEntry &b) {
        return (documentID == b.documentID) && (chromosome == b.chromosome);
    }
};

struct MatchingPair {
    unsigned int referencePartIndex;
    unsigned int readPartIndex;
    unsigned char referenceChromosome;
    MatchingPair(unsigned int refPart, unsigned int readPart, unsigned char refChro) {
        referencePartIndex = refPart;
        readPartIndex = readPart;
        referenceChromosome = refChro;
    }

    bool operator==(const MatchingPair &b) const{
        return referencePartIndex == b.referencePartIndex &&
                readPartIndex == b.readPartIndex &&
                referenceChromosome == b.referenceChromosome;
    }
};

namespace std
{
    template<>
    struct hash<MatchingPair>
    {
        typedef MatchingPair argument_type;
        typedef std::size_t value_type;

        value_type operator()(argument_type const& s) const
        {
            value_type const h1 ( std::hash<unsigned int>()(s.referencePartIndex) );
            value_type const h2 ( std::hash<unsigned int>()(s.readPartIndex) );
            value_type const h3 ( std::hash<unsigned int>()(s.referenceChromosome) );
            return (h1 ^ (h2 << 1)) ^ (h3 << 1);
        }
    };
}

class DocumentSet {
private:
    std::hash<unsigned int> hashFunction;
    
public:
    // Probably this way of separating the data structures instead of putting it into one class
    // is more cache-efficient, depending on the compiler

    // Saving the real string uses too much space, so we have to save the hash values of the
    // q-Grams
    std::vector<std::vector<unsigned int>> referenceSignatures, readSignatures;
    std::vector<unsigned int> referenceDocumentPositions, readDocumentPositions;
    std::vector<seqan::Dna5String> originalReadDocuments;
    std::vector<int> hashFunctions;
    std::unordered_map<unsigned int, std::vector<HashBucketEntry>> hashTable;

    std::unordered_set<MatchingPair> matchingPairs;
    std::vector<unsigned int> unmatchedReads;

    // This is a hash operating on a 32-bit string
    // Each char can so represent the unsigned int (also 32 bit) row hash value
    std::hash<std::u32string> bandHashFunction;

    size_t qGramSize;
    size_t signatureLength; // In other words: hash function count
    size_t bandSize; // For Locality-Sensitive Hashing
    size_t packageSize; // Number of reads in one package, that are checked for collisions and then send to the semi-global aligner 
    size_t limit; // The number of combinations that is at most added
    size_t limit_skip; // If the number of combinations in a q-gram is bigger than this value, no combination is added
    std::mt19937 rngLimit; // random number genetraotr
    QGramType type;

    void createHashFunctions();
    bool loadHashFunctions(std::istream &is);
    bool saveHashFunctions(std::ostream &os) const;
    bool loadReferenceSignatures(std::istream &is);
    bool saveReferenceSignatures(std::ostream &os) const;
    void addDocument(const seqan::Segment<const seqan::Dna5String> &document, unsigned int documentID,
                     std::vector< std::vector<unsigned int> > &documentSignaure);
    template <class IupacStringContainer>
    void addDocument(const IupacStringContainer &chromosome, size_t beginIndex, size_t windowSize,
                     unsigned int documentID, std::vector< std::vector<unsigned int> > &documentSignaure, bool useReverseComplement);

    template <class QGramContainer>
    void createSignatures(unsigned int documentID, QGramContainer &documentQGramSet,
                                       std::vector< std::vector<unsigned int> > &documentSignature);
    void findSimilarDocuments();
    void fillHashMapWithReferenceParts(unsigned int chromosome);
    void readFastaReferenceGenome(const seqan::Dna5String &genome, size_t windowSize, unsigned int chromosome);
    template <class MatchedFunction, class UnmatchedFunction>
        void readFastQReads(seqan::SequenceStream &reads, Timer &timer, MatchedFunction operationMatched, UnmatchedFunction operationUnmatched);
    template<class IupacStringContainer>
        void readReference(unsigned int chromosomeID, const IupacStringContainer &chromosomeString, size_t length, size_t windowSize);
    template <class Iterator> // an iterator pointing to seqan::Dna or seqan::Dna5
        static unsigned qGramHashFucntion(Iterator begin, Iterator end);

    void readFastaReferenceGenomeWithOffset(const seqan::Dna5String &genome, size_t windowSize, int &refNumber,
                                            size_t offset);


    void saveMatchingPairs(std::ostream &os) const;
    void loadMatchingPairs(std::istream &is);
    const std::unordered_set<MatchingPair> &getMatchingPairs() const {
        return matchingPairs;
    }

    DocumentSet();
    ~DocumentSet();
};


// =====================================================================================
// ===== I M P L E M E N T A T I O N S
// =====================================================================================

/**
 *  @param chromosome a object that overloads the operator[] and returns a iupac base.
 */
template <class IupacStringContainer>
void DocumentSet::addDocument(const IupacStringContainer &chromosome, size_t beginIndex, size_t windowSize,
                              unsigned int documentID, std::vector< std::vector<unsigned int> > &documentSignature, bool useReverseComplement) {
    // number of variants for a iupac symbol
    const static auto nv = [](unsigned x) {return (x & 1)
                                                + ((x & 2) >> 1)
                                                + ((x & 4) >> 2)
                                                + ((x & 8) >> 3);};
    const static unsigned long long NUMBER_OF_VARIANTS[16] = {nv(0), nv(1), nv(2), nv(3), nv(4), nv(5), nv(6), nv(7),
                                                              nv(8), nv(9), nv(10),nv(11),nv(12),nv(13),nv(14),nv(15)};


    std::unordered_set<unsigned int> documentQGramSet(3 * qGramSize); // initial size to avoid collions
    size_t endPos = beginIndex + windowSize - qGramSize + 1;



    // allocate space for findCombinationsOfQGram()
    std::vector<seqan::Dna> curQGram(qGramSize); // current q-gram
    std::vector<char> curVariant(qGramSize+1); // current variant stack. The +1 is important to avoid seg fault...
    typedef typename std::vector<seqan::Dna>::iterator OpIt;

    unsigned long long combinations = 1; // the current number of combinations
    for (size_t i = beginIndex; i < beginIndex + qGramSize; i++) {
        combinations *= NUMBER_OF_VARIANTS[(size_t) chromosome[i]];
    }

    // the normal operation to add all combinations
    auto operation = [&](OpIt beginIt, OpIt endIt) {
        documentQGramSet.insert(qGramHashFucntion(beginIt, endIt));
    };

    // the operation that adds on average only "limit" combinations
    const double rngMax = rngLimit.max();
    const double limitAsFloat = limit;
    auto operation_limit = [&documentQGramSet, this, rngMax, &combinations, limitAsFloat](OpIt beginIt, OpIt endIt) {
        double r = this->rngLimit();
        if (r / rngMax < limitAsFloat / combinations) {
            documentQGramSet.insert(qGramHashFucntion(beginIt, endIt));
        }
    };

    for (size_t position = beginIndex; position < endPos; position++) {
        if (combinations <= limit) {
            // add all combinations
            findCombinationsOfQGram(chromosome, position, qGramSize, curQGram, curVariant, operation);
        } else if (combinations <= limit_skip) {
            // add a random set of (on average) "limit" combinations
            findCombinationsOfQGram(chromosome, position, qGramSize, curQGram, curVariant, operation_limit);
        } else {
            // skip this q-gram (because there are too many combinations)
        }

        // update the number of combinations for the next q-gram
        if (position + 1 < endPos) { // true, if this is not the last iteration
            combinations /= NUMBER_OF_VARIANTS[(size_t) chromosome[position]]; // there is no remainder for this division!
            combinations *= NUMBER_OF_VARIANTS[(size_t) chromosome[position + qGramSize]];
        }
    }

    createSignatures(documentID, documentQGramSet, documentSignature);
}


/** Find all variant combinations of the substring "iupacString[beginPos...beginPos+qGramSize]"
 *  and execute the function operation on the founds.
 *  @param curQGram arbitrary space of length qGramSize
 *  @param curVariant arbitrary space of length qGramSize+1
 *  @param operation a object that overloads operator()(Iterator, Iterator) with Iterator = std::vector<seqan::DNA>::iterator
 */
template <class IupacString, class Function>
void findCombinationsOfQGram(IupacString &iupacString, size_t beginPos, size_t qGramSize,
                             std::vector<seqan::Dna> &curQGram, std::vector<char> &curVariant, Function& operation) {
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
        seqan::Iupac iupacSymbol = iupacString[beginPos + curIndex];
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

template <class Iterator> // an iterator pointing to seqan::Dna or seqan::Dna5
unsigned DocumentSet::qGramHashFucntion(Iterator begin, Iterator end) {
    static_assert(std::is_same<typename std::remove_reference<decltype(*begin)>::type, seqan::Dna>::value ||
                  std::is_same<typename std::remove_reference<decltype(*begin)>::type, seqan::Dna5>::value,
                  "Iterator is pointing to the wrong type");
    unsigned result = 0;
    for (; begin != end; ++begin) {
        result = result * 5 + unsigned(*begin);
    }
    return result;
}

template <class QGramContainer>
void DocumentSet::createSignatures(unsigned int documentID, QGramContainer &documentQGramSet,
                                   std::vector< std::vector<unsigned int> > &documentSignature) {
    std::vector<unsigned int> &signature = documentSignature.at(documentID);

    // Create signatureLength elements initialized with the maximum unsigned integer value
    signature.resize(signatureLength, std::numeric_limits<unsigned int>::max());

    for (unsigned int qgram : documentQGramSet) {
        unsigned int qGramHash = hashFunction(qgram);
        for (size_t hash = 0, end = signatureLength; hash < end; hash++) {
            // Create a normal hash value and XOR it with one of the randomly
            // generated integers. This will serve well enough as
            // "another" hash function
            unsigned int hashValue = qGramHash ^ hashFunctions[hash];
            signature[hash] = hashValue < signature[hash] ? hashValue : signature[hash];
        }
    }
}

template <class IupacStringContainer>
void DocumentSet::readReference(unsigned int chromosomeID,
                                const IupacStringContainer &chromosomeString,
                                size_t length,
                                size_t windowSize) {
    size_t refNumber;
    for (size_t offset: {(size_t) 0, windowSize / 2}) {
        for (size_t position = offset; position < length; position += windowSize, ++refNumber) {
            referenceDocumentPositions.push_back(position);
            referenceSignatures.push_back(std::vector<unsigned int>());
            addDocument(chromosomeString, position, windowSize, refNumber, referenceSignatures);
        }
    }

    fillHashMapWithReferenceParts(chromosomeID);
}

template <class MatchedFunction, class UnmatchedFunction>
void DocumentSet::readFastQReads(seqan::SequenceStream &reads, Timer &timer, MatchedFunction operationMatched, UnmatchedFunction operationUnmatched) {
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
//    seqan::readAll(ids, seqs, reads);
    //   size_t readNumber = seqan::length(seqs);

//    documents.resize(readNumber);
//    documentTypes.resize(documentTypes.size() + readNumber, READ);
    timer.printAndRestartTimer("Reading of reference and reads");

    
//    unsigned int packageCount = readNumber / packageSize;
//    for (size_t package = 0; package < packageCount; package++){
    size_t readCounter = 0;
    while (!seqan::atEnd(reads)) {
        seqan::clear(seqs);
        seqan::clear(ids);
        seqan::readBatch(ids, seqs, reads, packageSize);

        readSignatures.clear();
        readDocumentPositions.clear();
        originalReadDocuments.clear();

        readSignatures.resize(packageSize * 2);
        readDocumentPositions.resize(packageSize * 2);
        originalReadDocuments.resize(packageSize * 2);

        size_t currentPackageSize = seqan::length(seqs);
//        unsigned int packageEndId = ((package == (packageCount - 1))? readNumber : ((package + 1) * packageSize));

 //       std::cout << packageEndId << std::endl;
//#pragma omp parallel for
        size_t i = 0;
        for (; i < currentPackageSize; i++) {
            readDocumentPositions[i] = readCounter;
            originalReadDocuments[i] = seqs[i];
            //std::cout << originalReadDocuments.at(i) << std::endl;
            addDocument(seqs[i], i, readSignatures);
//#pragma omp atomic
            readCounter++;
        }
        seqan::reverseComplement(seqs);
        for (; i < currentPackageSize * 2; i++) {

            readDocumentPositions[i] = readCounter;
            originalReadDocuments[i] = seqs[i - currentPackageSize];
            //std::cout << originalReadDocuments.at(i) << std::endl;
            addDocument(seqs[i - currentPackageSize], i, readSignatures);
//#pragma omp atomic
            readCounter++;
        }
        findSimilarDocuments();

        operationMatched(matchingPairs);
        operationUnmatched(unmatchedReads);

        std::cout << "Count matching pairs: " << matchingPairs.size() << std::endl;
        std::cout << "Count not matching pairs: " << unmatchedReads.size() << std::endl;
        matchingPairs.clear();
        unmatchedReads.clear();
    }

//    std::cout << "Reads: " << readNumber << std::endl;
}


#endif // PG583_DOCUMENTSET_H
