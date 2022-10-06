//
// Created by daniel on 05/10/22.
//

#ifndef ABMCMC_SAMPLINGPROBLEM_H
#define ABMCMC_SAMPLINGPROBLEM_H

#include <functional>

template<class DISTRIBUTION, class SAMPLEUSER, class RESULT = std::result_of_t<SAMPLEUSER(std::function<typename DISTRIBUTION::domain_type()>)>>
class SamplingProblem {
public:
    typedef RESULT result_type;

    SAMPLEUSER user; // takes a producer of samples (from the distribution) and returns a result
    DISTRIBUTION distribution;


};


#endif //ABMCMC_SAMPLINGPROBLEM_H
