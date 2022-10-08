//
// Created by daniel on 07/10/22.
//

#ifndef ABMCMC_ASYNCVECTOR_H
#define ABMCMC_ASYNCVECTOR_H

#include <vector>
#include <future>

template<class PRODUCER, class OUTPUT = std::result_of_t<PRODUCER()>>
std::vector<OUTPUT> generateVectorAsync(int size, PRODUCER producerFunction) {
    if(size == 1) return {producerFunction()};
    std::vector<std::future<OUTPUT>> futures;
    for(int thread = 0; thread < size; ++thread) {
        futures.push_back(std::async(std::launch::async, producerFunction));
    }
    std::vector<OUTPUT> results;
    for(int thread = 0; thread < size; ++thread) {
        results.push_back(futures[thread].get());
    }
    return results;
}

#endif //ABMCMC_ASYNCVECTOR_H
