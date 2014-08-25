#ifndef PG583_TIMER_H
#define PG583_TIMER_H

#include <chrono>
#include <iostream>

// Example for usage:
//     Timer<> timer;
//     /* do something */
//     std::cout << timer << std::endl;
template <class Duration = std::chrono::milliseconds>
class Timer {
public:
    typedef std::chrono::high_resolution_clock::time_point TimePoint;
private:
    TimePoint timeStart;
    Duration duration;
public:
    Timer() :
        duration(Duration::zero()) {
        startTimer();
    }

    TimePoint now() const {
        return std::chrono::high_resolution_clock::now();
    }

    void startTimer() {
        timeStart = now();
    }

    Duration stopTimer() {
        TimePoint timeNow = now();
        duration += std::chrono::duration_cast<Duration>(timeNow - timeStart);
        timeStart = timeNow;
        return duration;
    }

    void resetTimer() {
        duration = Duration(Duration::zero());
        startTimer();
    }

    Duration getDuration() const {
        return duration;
    }

};

template <class Duration>
std::ostream& operator<<(std::ostream &o, Timer<Duration> &timer) {
    o << timer.stopTimer().count();
    return o;
}

#endif // PG583_TIMER_H
