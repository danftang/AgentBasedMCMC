//
// Created by daniel on 20/10/2021.
//

#ifndef GLPKTEST_COMPOUNDLOGGER_H
#define GLPKTEST_COMPOUNDLOGGER_H

#include <tuple>

template<typename...LOGGERS>
class CompoundLogger {
    std::tuple<LOGGERS...>  loggers;

    // Loggers should accept the operator << on samples
    CompoundLogger(LOGGERS &&... logrs): loggers(logrs...) {

    }

    template<class SAMPLE>
    void operator <<(const SAMPLE &sample) {
        unpackLoggers(sample, std::index_sequence_for<LOGGERS...>());
    }

    template<class SAMPLE, size_t... Indices>
    void unpackLoggers(const SAMPLE &sample, std::index_sequence<Indices...> indx) {
        ((std::get<Indices>(loggers) << sample),...);
    }

};


#endif //GLPKTEST_COMPOUNDLOGGER_H
