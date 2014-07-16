#ifndef TEST_H
#define TEST_H

#include <iostream>

class Test
{
public:
    Test();

    bool runTest();

    bool test1();

    bool test2();

    bool test3();

    bool test4();

    static bool testJQ();

    bool denseTest(std::string id, std::string reference, std::string read, std::string cigar);
};

#endif // TEST_H
