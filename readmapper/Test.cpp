/*#include "Test.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>

#include "align/Variant.h"
#include "align/SemiGlobalAligner.h"
#include <seqan/sequence.h>
#include <seqan/file.h>

Test::Test() {
}

bool Test::denseTest(std::string id, std::string reference, std::string read, std::string cigar) {

    int maxErrors = 5;
    std::vector<Chromosome> chromosomes;

    std::cout << "Test " << id;
    std::cout << " testing "
              << reference
              << " wich read " << read
              << " and cigar " << cigar << std::endl;

    // Generate Reference Object
    Chromosome chromosome("testChromosome1", reference, std::vector<Variant>());

    // Initialize read patterns
    ReadString pattern(read);
    auto chromosomeLength = seqan::length(chromosome.data);
    auto startPos = 1;
    auto endPos = chromosomeLength;

    // Voodoo
    VariantIndex variantIndex(chromosome.variants, chromosomeLength - 1);
    SemiGlobalAligner aligner(variantIndex);

    // Actual Alignment
    auto alignment = aligner.align(chromosome, startPos, endPos, pattern, maxErrors);
    std::string resultCigar = alignment.getCigar();

    std::cout << "awaited cigar " << cigar << std::endl;
    std::cout << "result cigar " <<  resultCigar << std::endl;

    // Test for equality
    if(resultCigar == cigar) {
        std::cout << "Test " << id << " OK!" << std::endl;

        return true;
    }

    std::cout << "Test " << id << " Fail!" << std::endl;
    return false;
}

bool Test::test1() {
    // Reference    "CCATACT GAACTG A CTAAC"
    // Read         "    ACTAGAA TG G CT"
    //                      I   D   X
    // Cigar        "3M1I3M1D5M"
    // Cigar ext    "3M1I3M1D2M1X2M"
    // Convertiert: "3M1D3M1I2M1X2M"

    return denseTest("1", "CCATACTGAACTGACTAAC", "ACTAGAATGGCT", "3M1D3M1I2M1X2M");
}

bool Test::test2() {
    // Reference "CACGATCA NN GACCGATACGTCCGA"
    // Read      "  CGATCA GA GACCGATA"
    // Cigar        "16M"

    return denseTest("2", "CACGATCANNGACCGATACGTCCGA", "CGATCAGAGACCGATA", "16M");
}

bool Test::test3() {
    // Reference    "AAAGGGAACCAAA"
    // Read         "       A CAA"
    // Cigar        "1M1D3M"
    // Convertiert  "1M1I3M"

    //              "2M1X2M


    return denseTest("3", "AAAGGGAACCAAA", "ACAA", "1M1I3M");
}

bool Test::test4() {
    // Reference "CACGATCA  GACCGATACGTCCGA"
    // Read      "  CGATCAGAGACCGATA"
    //              6M    2I8M
    // Cigar        "6M2I8M"
    // Convertiert: "6M2D8M"

    return denseTest("4", "CACGATCAGACCGATACGTCCGA", "CGATCAGAGACCGATA", "6M2D8M");
}

bool Test::runTest() {

    // Reference    CACGATCANNGACCGATACGTCCGA

    // Read         ATCANAGACCGATAC
    // Cigar        4M1P1I9M

    // Read         GATCANNGACCG
    // Cigar        R5M2P5M


    // Reference    AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG

    // Read         gggGTGTAACC-GACTAGgggg
    // Cigar        3S8M1D6M4S

    // Reference    AGCTAGCATCGTGTCGCCCGTCTAGCATACGCATGATCGACTGTCAGCTAGTCAGACTAGTCGATCGATGTG

    // Read         GTGTAACCC................................TCAGAATA (Intron)
    // Cigar        9M32N8M

    if (test1() & test2() & test3() & test4()) {
        std::cout << "All tests ran clean!" << std::endl;
        return true;
    }

    return false;
}

bool Test::testJQ() {
    ReferenceString ref  = "ACCGTATTCMRTCCCC";
    ReadString      read = "CGATACCTAAGTCC";

    Chromosome chromosome("0", ref, std::vector<Variant>());
    //    chromosome.variants.emplace_back((size_t)5, "ACGT");
    //    chromosome.variants.emplace_back((size_t)9, 3);
    //chromosome.variants.emplace_back((size_t)5, "ACGT");

    // Add some Variants
    chromosome.addInsertion((size_t) 5, "ACGT");            // Add Insertion
    chromosome.addDeletion((size_t) 5, (size_t) 3);         // Add Deletion Variant

    // Voodoo
    VariantIndex variantIndex(chromosome.variants, seqan::length(chromosome.data) - 1);
    SemiGlobalAligner aligner(variantIndex);

    // Alignemnt
    Alignment alignment = aligner.align(chromosome, 0, seqan::length(chromosome.data) - 1, read, 9999);

    std::cout << "Position: "<< alignment.getPosition() << std::endl;
    std::cout << "Variants: ";
    auto usedVariants = alignment.getUsedVariants();
    for (const Variant& v: usedVariants) {
        std::cout << v.getType() << ":" << v.getPosition() << " | ";
    }
    std::cout << std::endl;

    std::cout << "Cigar: " << alignment.getCigar() << std::endl;
    std::cout << "Align: " << alignment.getAlignmentString() << std::endl;

    return true;
}
*/
