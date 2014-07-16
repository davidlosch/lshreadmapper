#ifndef PG583_LOGGER_H
#define PG583_LOGGER_H

#include <iostream>

class Logger {
    std::ostream &infoStream;
    std::ostream &warnStream;
    std::ostream &errStream;
public:
    Logger(std::ostream &os = std::cerr);
    Logger(std::ostream &infoStream, std::ostream &warnStream, std::ostream &errStream);
    std::ostream& info(const char *prefix = "INFO: ");
    std::ostream& warn(const char *prefix = "WARNING: ");
    std::ostream& err(const char *prefix = "ERROR: ");
};

extern Logger logger;

#endif // PG583_LOGGER_H
