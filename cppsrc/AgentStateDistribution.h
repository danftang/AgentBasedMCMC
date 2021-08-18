//
// Created by daniel on 16/08/2021.
//

#ifndef GLPKTEST_AGENTSTATEDISTRIBUTION_H
#define GLPKTEST_AGENTSTATEDISTRIBUTION_H

// Represents a distribution over a given agent state
template<typename AGENT>
class AgentStateDistribution: public ConvexPMF {
public:
    AgentStateDistribution(const State<AGENT> &state, const ConvexPMF &oneDimensionalPMF);
};


#endif //GLPKTEST_AGENTSTATEDISTRIBUTION_H
