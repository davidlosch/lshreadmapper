#ifndef IUPACSTRING_H
#define IUPACSTRING_H

#include "seqan/sequence.h"
#include "seqan/seq_io.h"
#include <deque>

class IupacString
{
private:
    std::deque<char> data;

    unsigned char dna5ToIupacChar(int dna5);
    unsigned char asciiToIupac(char ascii);
public:
    IupacString();
    IupacString(seqan::String<seqan::Dna5String> &string);

    size_t size() const;
    seqan::Iupac operator[](size_t) const;
    void addSNP(size_t pos, seqan::Dna base);
};

#endif // IUPACSTRING_H
