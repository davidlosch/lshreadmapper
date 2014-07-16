#include "align/VariantIndex.h"
#include <vector>

using std::vector;

VariantIndex::VariantIndex(vector<Variant> &variantList, size_t highestPosition) :
    variantList(variantList) {
    // initialize private members
    size_t blockCount = (highestPosition >> 6) + 1; // one block for each 64 positions in reference
    existsVariant = vector<unsigned long long>(blockCount);
    listPositions = vector<size_t>(blockCount);
    listOffset = vector<vector<size_t> >(blockCount, vector<size_t>(0));

    // fill existence bitvector with zeroes
    for (size_t i = 0; i < blockCount; i++) {
        existsVariant[i] = 0;
    }

//    // sort variants by their position
//    std::sort(variantList.begin(), variantList.end(), [](Variant a, Variant b){
//        return a.getPosition() < b.getPosition();
//    });

    // iterate over variant list
    size_t oldBlock = 0;
    size_t oldPosition = 0;
    size_t oldOffset = 0;
    size_t variantsInCurrentBlock = 0;

    for (size_t i = 0; i < variantList.size(); i++) {
        size_t currentPosition = variantList[i].getPosition();
        assert(currentPosition >= oldPosition);
        size_t currentBlock = currentPosition >> 6;
        size_t currentOffset = currentPosition & 0x0000003f;
        if (currentBlock > oldBlock) {
            // shift bitvector for remaining number of positions in this block
            existsVariant[oldBlock] <<= 63 - oldOffset;
            oldOffset = 0;
            // new block entered: reset variant counter for current block and set start indes
            listPositions[currentBlock] = i;
            variantsInCurrentBlock = 0;
            oldBlock = currentBlock;
        }

        if (currentPosition > oldPosition) {
            // shift bitvector by difference of old and new offset
            existsVariant[currentBlock] <<= (currentOffset - oldOffset);
            // mark current offset for existing variant
            existsVariant[currentBlock]++;
            oldOffset = currentOffset;
//            printf("check 1 :: ");
            // position changed: add a new entry into offset list
//            printf("listOffset.size() = %u, currentBlock = %i, variants=%u", listOffset.size(), currentBlock, variantsInCurrentBlock);
            listOffset[currentBlock].push_back(variantsInCurrentBlock);
//            printf(" :: check 2\n");
            oldPosition = currentPosition;
        }
        variantsInCurrentBlock++;

    }
    // bitvector of last block has to be shifted for missing number of positions in it
    existsVariant[oldBlock] <<= 63 - oldOffset;
    ;
}

bool VariantIndex::isVariantAtPosition(size_t position) const {
    size_t block = position >> 6;  // = position/64
    size_t offset = position & 0x0000003f; // = position%64
    return (existsVariant[block] >> (63 - offset)) & 0x00000001;
}

unsigned int VariantIndex::getVariantsAtPosition(size_t position) const {
    size_t block = position >> 6;  // = position/64
    size_t offset = position & 0x0000003f; // = position%64
    if ((existsVariant[block] >> (63 - offset)) & 0x00000001) {
        unsigned int start = listPositions[block] + listOffset[block][popcount(existsVariant[block] >> (64 - offset))];
        return start;
    } else {
        return 0;
    }
}

const Variant& VariantIndex::getVariant(size_t position) const {
    return variantList[position];
}

size_t VariantIndex::getVariantCount() const {
    return variantList.size();
}

unsigned int VariantIndex::popcount(unsigned long long x) const {
    // implementation taken from "http://en.wikipedia.org/wiki/Hamming_weight"
    x -= (x >> 1) & 0x5555555555555555;                             //put count of each 2 bits into those 2 bits
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333); //put count of each 4 bits into those 4 bits
    x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;                        //put count of each 8 bits into those 8 bits
    return (x * 0x0101010101010101)>>56;                            //returns left 8 bits of x + (x<<8) + (x<<16) + ...
}
