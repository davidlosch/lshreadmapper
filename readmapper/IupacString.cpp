#include "IupacString.h"
#include "seqan/sequence.h"

// Deprecated!!!
// Klasse nicht mehr verwenden. Seqan stellt eine Methode bereit, die bereits Packed-Strings beinhaltet!

// Verwendung durch seqan::String<seqan::Iupac, seqan::Packed<>> nameDesStrings = "ACGT";
IupacString::IupacString(){

}

    IupacString::IupacString(seqan::String<seqan::Dna5String> &string){
        // Leider eine extremst hässliche Lösung, aber funktioniert.
        std::vector<char> chars(seqan::length(string)); //char chars[seqan::length(string)];
        for (size_t i = 0; i < seqan::length(string); i = i + 2){
            seqan::assign(chars, string[i]);
            unsigned char firstBlock = IupacString::asciiToIupac(chars[0]);
            unsigned char secondBlock = 0;
            if (i+1 < seqan::length(string)){
                seqan::assign(chars, string[i+1]);
                secondBlock = IupacString::asciiToIupac(chars[0]);
            }
            firstBlock = firstBlock << 4;
            unsigned char combinedBase = (firstBlock | secondBlock);
            data.push_back(combinedBase);
        }
    }

    unsigned char IupacString::asciiToIupac(char ascii){
        unsigned char base;
        switch (ascii) {
        case 'A':
            base = 1 << 1;
            break;
        case 'C':
            base = 1 << 2;
            break;
        case 'T':
            base = 1;
            break;
        case 'G':
            base = 1 << 3;
            break;
        default:
            base = 0;
            break;
        }
        return base;
    }

    unsigned char IupacString::dna5ToIupacChar(int dna5){
        unsigned char base;
        switch (dna5) {
        case 0:
            std::cout << "Case 0" << std::endl;
            base = 1 << 1;
            break;
        case 1:
            std::cout << "Case 1" << std::endl;
            base = 1 << 2;
            break;
        case 2:
            std::cout << "Case 2" << std::endl;
            base = 1 << 3;
            break;
        case 3:
            std::cout << "Case 3" << std::endl;
            base = 1;
            break;
        default:
            base = 0;
            break;
        }
        return base;
    }


size_t IupacString::size() const{
    // Query, if in the last Element of the data structure only the first 4 bits are set
    if ((data.back() & 15) == 0){
        return data.size()*2 - 1;
    }
    else {
        return data.size()*2;
    }
}

seqan::Iupac IupacString::operator[](size_t pos) const{
    unsigned char currentBase;
    if (pos % 2 == 0){
        currentBase = (data[pos/2] & 240) >> 4;
    }
    else {
        currentBase = data[pos/2] & 15;
    }
    int baseInt = (int) currentBase;
    seqan::Iupac iupacChar(baseInt);
    return iupacChar;
}

void IupacString::addSNP(size_t pos, seqan::Dna base){
    unsigned char currentBlock = data[pos/2];

    // Get both parts of the 8-bit char by bitwise AND
    unsigned char firstPart = currentBlock & 240;
    unsigned char secondPart = currentBlock & 15;

    unsigned char iupacBase = dna5ToIupacChar(base.value);

    // If the given position is even, the first part has to be altered
    if (pos % 2 == 0){
        firstPart = (firstPart | ((iupacBase & 15) << 4));
    }
    // otherwise the second
    else{
        secondPart = (secondPart | (iupacBase & 15));
    }
    currentBlock = (firstPart | secondPart);
    data[pos/2] = currentBlock;
}

