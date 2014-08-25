#include "align/Alignment.h"
#include "align/Variant.h"

Alignment::Alignment() :
    cigar(),
    alignmentString(),
    usedVariants(),
    startVariant(Variant::none),
    startOffset(0),
    errorCount(0),
    cigarInvalid(false),
    reversedCigar(),
    reversedAlignmentString(),
    currentCharCounter(0),
    currentChar(0) {
}

std::string Alignment::getCigar() {
    if (cigarInvalid) {
        // flush current char
        flush();
        currentChar = 0;
        currentCharCounter = 0;
        // rebuild cigar, if it has been changed
        cigar = std::string(reversedCigar.rbegin(), reversedCigar.rend());
        cigarInvalid = false;
    }
    return cigar;
}

ReferenceString Alignment::getAlignmentString() {
//    return reversedAlignmentString;
    if (alignmentStringInvalid) {
        // rebuild alignment string, if it has been changed
        seqan::clear(alignmentString);
        int last = seqan::length(reversedAlignmentString) - 1;
        for (int i = 0; i <= last; i++) {
           seqan::appendValue(alignmentString,reversedAlignmentString[last-i]);
        }
        alignmentStringInvalid = false;
    }
    return alignmentString;
}

const std::set<std::reference_wrapper<const Variant>>& Alignment::getUsedVariants() {
    return usedVariants;
}

size_t Alignment::getPosition() {
    return position;
}

const Variant& Alignment::getStartVariant() {
    return startVariant;
}

int Alignment::getStartOffset() {
    return startOffset;
}

size_t Alignment::getErrorCount() {
    return errorCount;
}

void Alignment::addCigarChar(char c) {
    if (c == currentChar) {
        // if added char is equal to last inserted, then just increment counter
        currentCharCounter++;
    } else {
        // write back last inserted char and its number of occurences (in reversed order)
        flush();
        // register new char
        currentChar = c;
        currentCharCounter = 1;
    }
    cigarInvalid = true;
}

void Alignment::addAlignmentChar(ReferenceChar c) {
    seqan::appendValue(reversedAlignmentString,c);
    alignmentStringInvalid = true;
}

void Alignment::addUsedVariant(const Variant &var) {
    if (&var != &Variant::none) {
        usedVariants.insert(var);
    }
}

void Alignment::setPosition(size_t pos) {
    position = pos;
}

void Alignment::setStartVariant(const Variant &var) {
    startVariant = var;
}

void Alignment::setStartOffset(int off) {
    startOffset = off;
}

void Alignment::flush() {
    if(currentCharCounter > 0) {
        reversedCigar += currentChar;
        std::string count = std::to_string(currentCharCounter);
        reversedCigar += std::string(count.rbegin(), count.rend());
    }
}

void Alignment::incrementErrorCount() {
    errorCount++;
}
