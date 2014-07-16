#include "Logger.h"

Logger::Logger(std::ostream &os/* = std::cerr*/) :
    Logger(os, os, os) {}
Logger::Logger(std::ostream &infoStream, std::ostream &warnStream, std::ostream &errStream) :
    infoStream(infoStream),
    warnStream(warnStream),
    errStream(errStream)
{}
std::ostream& Logger::info(const char *prefix/* = "INFO: "*/) {
    return infoStream << prefix;
}
std::ostream& Logger::warn(const char *prefix/* = "WARNING: "*/) {
    return warnStream << prefix;
}
std::ostream& Logger::err(const char *prefix/* = "ERROR: "*/) {
    return errStream << prefix;
}

Logger logger;


