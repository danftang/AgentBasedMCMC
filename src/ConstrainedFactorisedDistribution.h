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
#include <cmath>
#include "EqualityConstraints.h"
#include "SparseFunction.h"
#include "ABM.h"

template<class COEFF> class TableauNormMinimiser;

template<typename DOMAIN, typename CONSTRAINTCOEFF = typename subscript_operator_traits<DOMAIN>::decay_type>
class ConstrainedFactorisedDistribution {
public:
    typedef CONSTRAINTCOEFF coefficient_type;
    typedef SparseFunction<std::pair<double,bool>,const DOMAIN &> function_type;
    typedef DOMAIN domain_type;

    std::vector<function_type>           factors;
    EqualityConstraints<CONSTRAINTCOEFF> constraints;       // linear constraints

    void factoriseConstraints();

    void addFactor(SparseFunction<std::pair<double,bool>,const DOMAIN &> factor) {
        factors.push_back(std::move(factor));
    }

    template<class INDEXDEPENDENCIES = std::initializer_list<int>>
    void addFactor(std::function<std::pair<double,bool>(const DOMAIN &)> func, const INDEXDEPENDENCIES &indexDependencies) {
        factors.emplace_back(std::move(func), indexDependencies);
    }

    void addConstraint(EqualityConstraint<CONSTRAINTCOEFF> constraint) {
        constraints.push_back(std::move(constraint));
    }

    void addConstraints(EqualityConstraints<CONSTRAINTCOEFF> constraints) {
        for(auto &constraint : constraints) {
            this->constraints.push_back(std::move(constraint));
        }
    }

    bool isFeasible(const DOMAIN &X) {
        if(!constraints.template isValidSolution(X)) return false;
        for(const auto &factor: this->factors) {
            if(factor(X).second == false) return false;
        }
        return true;
    }


    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &operator *=(const ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &other) {
        constraints += other.constraints;
        factors.reserve(factors.size() + other.factors.size());
        factors.insert(factors.end(), other.factors.begin(), other.factors.end());
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

    ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> operator *(ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF> &&factoredDist) && {
        factoredDist *= *this;
        return std::move(factoredDist);
    }


    double operator ()(const DOMAIN &X) const { return Pexact(X); }

    double logPexact(const DOMAIN &X) const {
        double logP = 0.0;
        if(!constraints.isValidSolution(X)) return -std::numeric_limits<double>::infinity();
        for(const auto &factor: factors) {
            std::pair<double,bool> factorVal = factor(X);
            if(factorVal.second == false) return -std::numeric_limits<double>::infinity();
            logP += factorVal.first;
        }
        return logP;
    }

    // use logPexact if the total probability may be too small for double.
      double Pexact(const DOMAIN &X) const { return std::exp(logPexact(X)); }

//    double logPwidened(const DOMAIN &X) const {
//        double logP = 0.0;
//        double distanceToValidHyperplane = 0.0;
//        for(const auto &constraint : constraints) {
//            distanceToValidHyperplane += fabs(constraint.coefficients * X - constraint.constant);
//        }
//        logP += -ABM::kappa * distanceToValidHyperplane;
//        for(const auto &factor: factors) {
//            logP += log(factor(X).first);
//        }
//        return logP;
//    }

    void sanityCheck(const DOMAIN &testVector) const {
        for(const auto &factor: factors) {
            factor.sanityCheck(testVector);
        }
    }

//    template<typename = std::is_same<ModelState<typename DOMAIN::agent_type>, DOMAIN>>
//    auto toTrajectoryDistribution(int time) const {
//        ConstrainedFactorisedDistribution<Trajectory<typename DOMAIN::agent_type>> trajDistribution;
//        trajDistribution.factors = trajectoryFactors(time);
//        trajDistribution.constraints = trajectoryConstraints(time);
//        return trajDistribution;
//    }

//    template<class NEWDOMAIN, typename = decltype(typename DOMAIN::template linearMapTo<NEWDOMAIN>())>
//    ConstrainedFactorisedDistribution<NEWDOMAIN> toDomain() {
//        ConstrainedFactorisedDistribution<NEWDOMAIN> trajDistribution;
//        trajDistribution.factors = trajectoryFactors(time);
//        trajDistribution.constraints = trajectoryConstraints(time);
//        return trajDistribution;
//
//    }

//    std::vector<SparseVec<CONSTRAINTCOEFF>> calculateBasis(int domainDimension) const {
//        std::vector<SparseVec<CONSTRAINTCOEFF>> basisVectors;
//        TableauNormMinimiser<CONSTRAINTCOEFF> constraintTableau(constraints);
//        std::vector<SparseVec<CONSTRAINTCOEFF>> uniqueBasisVectors = constraintTableau.getBasisVectors(domainDimension);
//        // double-up unique basis vectors for + and -
//        basisVectors.reserve(uniqueBasisVectors.size()*2);
//        for(int i=0; i < uniqueBasisVectors.size(); ++i) {
//            basisVectors.push_back(-uniqueBasisVectors[i]);
//            basisVectors.push_back(std::move(uniqueBasisVectors[i]));
//        }
//        return basisVectors;
//    }


protected:

    // Given a linear map from domain A->B (i.e. B = MA + C) we can transform a disrtibution over B into a distribution over A
    // ...special case is when M is binary and C=0, so B_i = sum_{j \in S(i)} A_j so we only need to supply S(i)

    // If this distribution's domain is ModelState<AGENT> then this
    // generates a set of factors on Trajectory<AGENT> domain, with the
    // constraints applied at the given time.
//    template<typename AGENT = typename DOMAIN::agent_type, typename = std::is_same<ModelState<AGENT>, DOMAIN>>
//    auto trajectoryFactors(int time) {
//        std::vector<SparseFunction<std::pair<double, bool>, const Trajectory<AGENT> &>> trajFactors;
//        for (const auto &modelStateFactor: factors) {
//            std::vector<int> trajectoryDependencies;
//            trajectoryDependencies.reserve(modelStateFactor.dependencies.size() * AGENT::actDomainSize);
//            for (int agentId: modelStateFactor.dependencies) {
//                for (int actId: State<AGENT>(time, agentId).forwardOccupationDependencies()) {
//                    trajectoryDependencies.push_back(actId);
//                }
//            }
//            trajFactors.emplace_back(
//                    [modelFactor = modelStateFactor, time](const Trajectory<AGENT> &trajectory) {
//                        return modelFactor(trajectory.temporaryPartialModelState(time,modelFactor.dependencies));
//                    },
//                    trajectoryDependencies
//            );
//        }
//        return trajFactors;
//    }


    // Transforms factors to a new domain subject to a linear transform from the new domain to the
    // old. A = MB
//    template<typename NEWDOMAIN>
//    auto transformedFactors() {
//        const std::vector<SparseVec<CONSTRAINTCOEFF>> &linearTransform = NEWDOMAIN::template linearDomainTransform<DOMAIN>();
//        std::vector<SparseFunction<std::pair<double, bool>, const NEWDOMAIN &>> newFactors;
//        for (const auto &oldFactor: factors) {
//            std::set<int> newDependencies;
//            for (int agentId: oldFactor.dependencies) {
//                for (int actId: ) {
//                    newDependencies.push_back(actId);
//                }
//            }
//            newFactors.emplace_back(
//                    [modelFactor = oldFactor, time](const Trajectory<AGENT> &trajectory) {
//                        return modelFactor(trajectory.temporaryPartialModelState(time,modelFactor.dependencies));
//                    },
//                    newDependencies
//            );
//        }
//        return newFactors;
//    }


//    template<typename AGENT = typename DOMAIN::agent_type, typename = std::is_same<ModelState<AGENT>,DOMAIN>>
//    EqualityConstraints<CONSTRAINTCOEFF> trajectoryConstraints(int time) {
//        EqualityConstraints<CONSTRAINTCOEFF> trajConstraints;
//        for(auto &modelStateConstraint: constraints) {
//            EqualityConstraint<CONSTRAINTCOEFF> &trajectoryConstraint = trajConstraints.emplace_back();
//            for(int nzi = 0; nzi < modelStateConstraint.coefficients.sparseSize(); ++nzi) {
//                auto actIds = State<AGENT>(time, modelStateConstraint.coefficients.indices[nzi]).forwardOccupationDependencies();
//                for(int actId: actIds) trajectoryConstraint.coefficients.insert(actId, modelStateConstraint.coefficients.values[nzi]);
//            }
//            trajectoryConstraint.constant = modelStateConstraint.constant;
//        }
//        return trajConstraints;
//    }
};


template<typename T>
std::ostream &operator <<(std::ostream &out, const ConstrainedFactorisedDistribution<T> &distribution) {
    for(int i=0; i<distribution.factors.size() && i<100; ++i) {
        out << "P(";
        auto deps = distribution.factors[i].dependencies;
        if(deps.size() > 0) {
            out << "X" << deps[0];
            for(int i=1; i<deps.size(); ++i) {
                out << ", X" << deps[i];
            }
        }
        out << ")";
    }
    if(distribution.factors.size() > 100) out << "... ...";
    out << "\nSubject to\n";
    out << distribution.constraints;
    return out;
}


template<class DOMAIN, class CONSTRAINTCOEFF>
void ConstrainedFactorisedDistribution<DOMAIN,CONSTRAINTCOEFF>::factoriseConstraints() {
    TableauNormMinimiser<CONSTRAINTCOEFF>::factoriseConstraints(*this);
}


#endif //ABMCMC_CONSTRAINEDFACTORISEDDISTRIBUTION_H
