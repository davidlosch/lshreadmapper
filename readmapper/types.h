#ifndef PG583_TYPES_H
#define PG583_TYPES_H

#include <seqan/sequence.h>

#define USE_PACKED_IUPAC_STRING
#ifdef USE_PACKED_IUPAC_STRING
typedef seqan::Iupac ReferenceChar;
typedef seqan::String<ReferenceChar, seqan::Packed<>> ReferenceString;
typedef seqan::Dna5/*Iupac*/ ReadChar;
typedef seqan::String<ReadChar/*, seqan::Packed<>*/> ReadString;
#else
#include <string>
typedef char ReferenceChar;
typedef char ReadChar;
typedef std::basic_string<ReferenceChar> ReferenceString;
typedef std::basic_string<ReadChar> ReadString;
#endif
typedef seqan::Value<ReferenceChar>::Type ReferenceCharValue;
typedef seqan::Value<ReadChar>::Type ReferenceCharValue;

#endif // PG583_TYPES_H
