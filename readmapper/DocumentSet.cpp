#include "DocumentSet.h"

#include <limits>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>

DocumentSet::DocumentSet() :
    rngLimit(time(NULL)) {
    type = Q_GRAM_CHAR;
    qGramSize = 10;
    signatureLength = 100;
    bandSize = 2;
    limit = 64;
    limit_skip = 256;
}

void DocumentSet::addDocument(const seqan::Segment<const seqan::Dna5String> &document, unsigned int documentID,
                              std::vector< std::vector<unsigned int> > &documentSignature) {
    // Just buffering the constant document length
   size_t docLength = seqan::length(document);
   std::vector<unsigned int> documentQGramSet;// = documents.at(documentID);

   // Working with single characters as Q-Gram elements is more straight forward..
   assert (type == Q_GRAM_CHAR);

//   currentDocument.reserve(docLength);
   for (size_t position = 0; position < docLength - qGramSize; position++) {
       size_t infixend = position + qGramSize;
       auto infix = seqan::infix(document, position, infixend);
       unsigned int hashValue = infix[0];
       for (size_t i = 1; i < qGramSize; i++) {
           hashValue = hashValue * 5 + (unsigned int) infix[i];
       }
       documentQGramSet.push_back(hashValue);
   }
   // Remove duplicates.
   std::sort(documentQGramSet.begin(), documentQGramSet.end());
   documentQGramSet.erase(std::unique(documentQGramSet.begin(), documentQGramSet.end()), documentQGramSet.end());
   createSignatures(documentID, documentQGramSet, documentSignature);
}



void DocumentSet::createHashFunctions() {
    std::mt19937 rng;
    rng.seed(0);
    hashFunctions.clear();
    for (size_t i = 0; i < signatureLength; i++) {
        hashFunctions.push_back(rng());
    }
}

bool DocumentSet::loadHashFunctions(std::istream &is)
{
    if (!is.good()) {
        return false;
    }

    is >> std::ws >> signatureLength;

    for (size_t i = 0; i < signatureLength; i++) {
        if (!is.good()) {
            return false;
        }

        unsigned int hash;
        is >> hash;
        hashFunctions.push_back(hash);
    }

    return true;
}

bool DocumentSet::saveHashFunctions(std::ostream &os) const {
    os << signatureLength << '\n';
    for (auto hash : hashFunctions) {
        os << hash << '\n';
    }
    return true;
}

//void DocumentSet::createSignatures(unsigned int documentID, std::vector<unsigned int> &documentQGramSet,
//                                   std::vector< std::vector<unsigned int> > &documentSignature) {
//    std::vector<unsigned int> &signature = documentSignature.at(documentID);

//    // Create signatureLength elements initialized with the maximum integer value
//    signature.resize(signatureLength, std::numeric_limits<int>::max());

//    for (unsigned int qgram : documentQGramSet) {
//        unsigned int qGramHash = hashFunction(qgram);
//        for (size_t hash = 0, end = signatureLength; hash < end; hash++) {
//            // Create a normal hash value and XOR it with one of the randomly
//            // generated integers. This will serve well enough as
//            // "another" hash function
//            unsigned int hashValue = qGramHash ^ hashFunctions[hash];
//            signature[hash] = hashValue < signature[hash] ? hashValue : signature[hash];
//        }
//    }
//}

void DocumentSet::findSimilarDocuments() {
    size_t docSize = readSignatures.size();
    
    unsigned int bandID = 0;

    for (unsigned int i = 0; i < docSize; i++) {
        std::vector<unsigned int> &signature = readSignatures.at(i);
        std::u32string band;
        bool foundAtLeastOneMatch = false;
       
        for (unsigned int &signatureValue : signature) {
            band.push_back((char32_t) signatureValue);

            if (band.length() == bandSize) {
                unsigned int bandHashValue = bandHashFunction(band);
                auto bucket = hashTable.find(bandHashValue);
                if (bucket != hashTable.end()) {
                    for (HashBucketEntry &entry : bucket->second) {
//                        if (entry.bandID == bandID){
//                            MatchingPair part = { entry.documentID, i, (unsigned char)entry.chromosome };
                            matchingPairs.emplace(entry.documentID, i, (unsigned char)entry.chromosome );
                            //std::cout << originalReadDocuments.at(i) << std::endl;
                            //std::cout << originalReferenceDocuments.at(entry) << std::endl << std::endl;
                            foundAtLeastOneMatch = true;
//                        }
                    }
                } 
                band.clear();
                bandID++;
            }
        }
        bandID = 0;
        if (!foundAtLeastOneMatch){
            unmatchedReads.push_back(i); 
        }
    }
    
    std::cout << matchingPairs.size() << " matching pairs" << std::endl;
}

void DocumentSet::saveMatchingPairs(std::ostream &os) const {
    unsigned int numMatches = matchingPairs.size();
    os.write(reinterpret_cast<char *>(&numMatches), sizeof(numMatches));

    // TODO
//    os.write(reinterpret_cast<const char *>(matchingPairs.data()), numMatches * sizeof(MatchingPair));
}

void DocumentSet::loadMatchingPairs(std::istream &is) {
    unsigned int numMatches;
    is.read(reinterpret_cast<char *>(&numMatches), sizeof(numMatches));
    // TODO
//    matchingPairs.resize(numMatches);
//    is.read(reinterpret_cast<char *>(matchingPairs.data()), numMatches * sizeof(MatchingPair));
}

bool DocumentSet::loadReferenceSignatures(std::istream &is) {
    if (!is.good()) {
        return false;
    }

    int numSignatures;
    is.read(reinterpret_cast<char *>(&signatureLength), sizeof(int));
    is.read(reinterpret_cast<char *>(&numSignatures), sizeof(int));

    for (int i = 0; i != numSignatures; ++i) {
        std::vector<unsigned> sig;
        sig.resize(signatureLength);
        is.read(reinterpret_cast<char *>(sig.data()), signatureLength * sizeof(unsigned));
        referenceSignatures.push_back(std::move(sig));
    }

    // TODO:
    // Save / Load chromosome ID to file
    fillHashMapWithReferenceParts(/* chromosome = */ 0);
    return is.good();
}

bool DocumentSet::saveReferenceSignatures(std::ostream &os) const {
    os.write(reinterpret_cast<const char *>(&signatureLength), sizeof(int));
    int numSignatures = referenceSignatures.size();
    os.write(reinterpret_cast<const char *>(&numSignatures), sizeof(int));

    for (const auto &sig : referenceSignatures) {
        os.write(reinterpret_cast<const char *>(sig.data()), sig.size() * sizeof(unsigned));
    }

    return os.good();
}

void DocumentSet::fillHashMapWithReferenceParts(unsigned int chromosome) {
    size_t docSize = referenceSignatures.size();
    unsigned int bandID = 0;

    for (unsigned int i = 0; i < docSize; i++) {
        std::vector<unsigned int> &signature = referenceSignatures.at(i);
        std::u32string band;
        for (unsigned int &signatureValue : signature) {
            band.push_back((char32_t) signatureValue);

            if (band.length() == bandSize) {
                unsigned int bandHashValue = bandHashFunction(band);
                std::vector<HashBucketEntry> &bucket = hashTable[bandHashValue];
                HashBucketEntry currentRefPart = {i, chromosome/*, bandID++*/};
                if (std::find(bucket.begin(), bucket.end(), currentRefPart) == bucket.end()) {

                    // Add new reference genome part to this bucket, if it was not found
                    bucket.push_back(currentRefPart);
                }
                band.clear();
            }
        }
        bandID = 0;
    }

//    originalReferenceDocuments.clear();
    referenceSignatures.clear();
//    referenceDocumentPositions.clear();
}

void DocumentSet::readFastaReferenceGenomeWithOffset(const seqan::Dna5String &genome, size_t windowSize,
                                                     int &refNumber, size_t offset) {
    for (size_t position = offset, end = seqan::length(genome);
         position < end;
         position += windowSize, ++refNumber) {
        size_t infixEnd = std::min(position + windowSize, end);
        auto window = seqan::infix(genome, position, infixEnd);
        referenceDocumentPositions.push_back(position);
//        originalReferenceDocuments.push_back(window);
//        referenceDocuments.push_back(std::vector<unsigned int>());
        referenceSignatures.push_back(std::vector<unsigned int>());
        addDocument(window, refNumber, referenceSignatures);
        //std::cout << "Adding part " << refNumber << std::endl;
    }

}

void DocumentSet::readFastaReferenceGenome(const seqan::Dna5String &genome, size_t windowSize, unsigned int chromosome ) {
    // Now create the reference genome parts as documents
    int refNumber = 0;
    readFastaReferenceGenomeWithOffset(genome, windowSize, refNumber, 0);
    readFastaReferenceGenomeWithOffset(genome, windowSize, refNumber, (size_t)(windowSize >> 1));
    fillHashMapWithReferenceParts(chromosome);
    std::cout << "Reference Parts: " << refNumber << std::endl;
}

DocumentSet::~DocumentSet() {

}
