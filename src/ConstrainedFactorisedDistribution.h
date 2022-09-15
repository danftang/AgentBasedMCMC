// A ConstrainedFactorisedDistribution represents a probability distribution over
// the template class T of the form
//
// P(X) = \prod_i P_i(X)
//
// subject to
//
// MX = C
// X \in {set of integers}
//
// for some set of sparse, widened functions P_i.
// The template class T should have a subscript operator.
//
// Created by daniel on 01/09/22.
//

#ifndef ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H
#define ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H

#include <vector>
#include <limits>
#include "EqualityConstraints.h"
#include "SparseFunction.h"

template<typename DOMAIN, typename CONSTRAINTCOEFF = typename subscript_operator_traits<DOMAIN>::base_type>
class ConstrainedFactorisedDistribution {
public:
    typedef CONSTRAINTCOEFF coefficient_type;
    typedef SparseFunction<std::pair<double,bool>,const DOMAIN &> function_type;
    typedef DOMAIN domain_type;

    std::vector<function_type>           logFactors;
    EqualityConstraints<CONSTRAINTCOEFF> constraints;        // linear constraints
//    int                                  domainDimension;


    void addFactor(SparseFunction<std::pair<double,bool>,const DOMAIN &> factor) {
        logFactors.push_back(std::move(factor));
    }

    void addConstraint(EqualityConstraint<CONSTRAINTCOEFF> constraint) {
//        this->domainDimension = std::max(this->domainDimension, constraint.maxCoefficientIndex());
        constraints.push_back(std::move(constraint));
    }

    double exactFactorValue(int factorIndex, const DOMAIN &X) const {
        std::pair<double,bool> factorVal = logFactors[factorIndex](X);
        return factorVal.second?factorVal.first:-std::numeric_limits<double>::infinity();
    }

    double widenedFactorValue(int factorIndex, const DOMAIN &X) const {
        return logFactors[factorIndex](X).first;
    }


    bool isFeasible(const DOMAIN &X) {
        if(!constraints.template isValidSolution(X)) return false;
        for(const auto &factor: this->logFactors) {
            if(factor(X).second == false) return false;
        }
        return true;
    }



    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &operator *=(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &other) {
        constraints += other.constraints;
        logFactors.reserve(logFactors.size() + other.logFactors.size());
        logFactors.insert(logFactors.end(), other.logFactors.begin(), other.logFactors.end());
//        if(other.domainDimension > domainDimension) domainDimension = other.domainDimension;
        return *this;
    }


    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &factoredDist) && {
        (*this) *= factoredDist;
        return std::move(*this);
    }

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &factoredDist) const & {
        ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> copyOfThis(*this);
        copyOfThis *= factoredDist;
        return copyOfThis;
    }

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &&factoredDist) const & {
        factoredDist *= *this;
        return std::move(factoredDist);
    }

    double logPexact(const DOMAIN &X) const {
        double logP = 0.0;
        if(!constraints.isValidSolution(X)) return -std::numeric_limits<double>::infinity();
        for(int i=0; i < this->logFactors.size(); ++i) {
            logP += this->exactFactorValue(i,X);
        }
        return logP;
    }

    double logPwidened(const DOMAIN &X) const {
        double logP = 0.0;
        double distanceToValidHyperplane = 0.0;
        for(const auto &constraint : constraints) {
            distanceToValidHyperplane += fabs(constraint.coefficients * X - constraint.constant);
        }
        logP += -ABM::kappa * distanceToValidHyperplane;
        for(int i=0; i < this->logFactors.size(); ++i) {
            logP += this->widenedFactorValue(i,X);
        }
        return logP;
    }

    void sanityCheck(const DOMAIN &testVector) const {
        for(const auto &factor: logFactors) {
            factor.sanityCheck(testVector);
        }
    }

    template<typename = std::is_same<ModelState<typename DOMAIN::agent_type>, DOMAIN>>
    auto toTrajectoryDistribution(int time) {
        ConstrainedFactorisedDistribution<Trajectory<typename DOMAIN::agent_type>> trajDistribution;
        trajDistribution.logFactors = trajectoryLogFactors(time);
        trajDistribution.constraints = trajectoryConstraints(time);
        return trajDistribution;
    }


protected:

    template<typename AGENT = typename DOMAIN::agent_type, typename = std::is_same<ModelState<AGENT>, DOMAIN>>
    auto trajectoryLogFactors(int time) {
        std::vector<SparseFunction<std::pair<double, bool>, const Trajectory<AGENT> &>> trajFactors;
        for (const auto &modelStateFactor: logFactors) {
            std::vector<int> trajectoryDependencies;
            trajectoryDependencies.reserve(modelStateFactor.dependencies.size() * AGENT::actDomainSize());
            for (int agentId: modelStateFactor.dependencies) {
                for (int actId: State<AGENT>(time, agentId).forwardOccupationDependencies()) {
                    trajectoryDependencies.push_back(actId);
                }
            }
            trajFactors.emplace_back(
                    [modelFactor = modelStateFactor, time](const Trajectory<AGENT> &trajectory) { // TODO: ensure that modelFactor is caprured by value
                        return modelFactor(trajectory.temporaryPartialModelState(time,modelFactor.dependencies));
                    },
                    trajectoryDependencies
            );
        }
        return trajFactors;
    }


    template<typename AGENT = typename DOMAIN::agent_type, typename = std::is_same<ModelState<AGENT>,DOMAIN>>
    EqualityConstraints<CONSTRAINTCOEFF> trajectoryConstraints(int time) {
        EqualityConstraints<CONSTRAINTCOEFF> trajConstraints;
        for(auto &modelStateConstraint: constraints) {
            EqualityConstraint<CONSTRAINTCOEFF> &trajectoryConstraint = trajConstraints.emplace_back();
            for(int nzi = 0; nzi < modelStateConstraint.coefficients.sparseSize(); ++nzi) {
                auto actIds = State<AGENT>(time, modelStateConstraint.coefficients.indices[nzi]).forwardOccupationDependencies();
                for(int actId: actIds) trajectoryConstraint.coefficients.insert(actId, modelStateConstraint.coefficients.values[nzi]);
            }
            trajectoryConstraint.constant = modelStateConstraint.constant;
        }
        return trajConstraints;
    }
};


template<typename T>
std::ostream &operator <<(std::ostream &out, const ConstrainedFactorisedDistribution<T> &distribution) {
    for(const auto &factor: distribution.logFactors) {
        out << "P(";
        if(factor.dependencies.size() > 0) {
            out << "X" << factor.dependencies[0];
            for(int i=1; i<factor.dependencies.size(); ++i) {
                out << ", X" << factor.dependencies[i];
            }
        }
        out << ")";
    }
    out << "\nSubject to\n";
    out << distribution.constraints;
    return out;
}

#endif //ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H
