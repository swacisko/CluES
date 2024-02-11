//
// Created by sylwester on 3/9/21.
//

#ifndef ALGORITHMSPROJECT_GLOBAL_H
#define ALGORITHMSPROJECT_GLOBAL_H

#include <csignal>
#include <mutex>
#include <sys/resource.h>
#include "Makros.h"

namespace Global{

    extern map<string,int64_t> counters;

    extern int max_runtime_in_seconds;

    /**
     * If true, then no logs will be written to [clog]
     */
    extern const bool disable_all_logs;

    /**
     *
     * @return true if TLE, false otherwise
     */
    extern bool checkTle();

    extern void startAlg();

    /**
     * @return number of seconds from the start of the algorithm
     */
    extern int secondsFromStart();


    /**
     * Sets stack size
     */
    void increaseStack();


    /**
     * If true, then
     */
    extern const bool CONTEST_MODE;
}

#endif //ALGORITHMSPROJECT_GLOBAL_H
