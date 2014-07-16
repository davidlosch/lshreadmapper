#ifndef PG583_TIMER_H
#define PG583_TIMER_H

#include <chrono>
#include <string>
#include <iostream>

class Timer {
private:
    std::chrono::high_resolution_clock::time_point timeStart, timeEnd;
    std::chrono::milliseconds duration;
public:
    Timer() {}

    inline void startTimer() {
        timeStart = std::chrono::high_resolution_clock::now();
    }

    inline void stopTimer() {
        timeEnd = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeStart);
    }

    inline void restartTimer() {
        stopTimer();
        timeStart = std::chrono::high_resolution_clock::now();
    }

    void printDuration(std::string str) {
        std::cout << str << ": " << std::fixed << duration.count() << std::endl;
    }

    inline void printAndRestartTimer(std::string str) {
        stopTimer();
        printDuration(str);
        startTimer();
    }

};


// a better timer class
// Example for usage:
//     JQTimer<> timer;
//     /* do something */
//     std::cout << timer << std::endl;
template <class Measurement = std::chrono::milliseconds>
class JQTimer {
public:
    typedef std::chrono::high_resolution_clock::time_point TimePoint;
private:
    TimePoint timeStart;
    Measurement duration;
public:
    JQTimer() :
        duration(0) {
        startTimer();
    }

    TimePoint now() const {
        return std::chrono::high_resolution_clock::now();
    }

    void startTimer() {
        timeStart = now();
    }

    void stopTimer() {
        TimePoint timeNow = now();
        duration += std::chrono::duration_cast<Measurement>(timeNow - timeStart);
        timeStart = timeNow;
    }

    void resetTimer() {
        duration = Measurement(0);
        startTimer();
    }

    Measurement getDuration() {
        stopTimer();
        return duration;
    }

};

template <class Measurement>
std::ostream& operator<<(std::ostream& o, JQTimer<Measurement>& timer) {
    o << timer.getDuration().count();
    return o;
}

#endif // PG583_TIMER_H
