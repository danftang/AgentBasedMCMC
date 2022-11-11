//
// Created by daniel on 04/10/22.
//

#ifndef ABMCMC_ABMPOSTERIOR_H
#define ABMCMC_ABMPOSTERIOR_H

#include "State.h"
#include "ConstrainedFactorisedDistribution.h"
#include "ABMLikelihood.h"
#include "ABMPrior.h"
#include "Basis.h"

template<class TRAJECTORY, class STARTSTATE>
class ABMPosterior: public ConstrainedFactorisedDistribution<TRAJECTORY> {
public:
    typedef typename TRAJECTORY::agent_type agent_type;

    ABMPrior<TRAJECTORY, STARTSTATE>    prior;
    ABMLikelihood<TRAJECTORY>           likelihood;
    Basis<TRAJECTORY>                   basis;

    ABMPosterior()=default;

    ABMPosterior(ABMPrior<TRAJECTORY,STARTSTATE> prior, ABMLikelihood<TRAJECTORY> likelihood):
            prior(std::move(prior)),
            likelihood(std::move(likelihood)) {
        init();
    }


    ABMPosterior(STARTSTATE startState, ABMLikelihood<TRAJECTORY> likelihood):
    prior(std::move(startState)),
    likelihood(std::move(likelihood)) {
        init();
    }

//    ABMPosterior(STARTSTATE startState, double pMakeObservation, double pObserveIfPresent, double kappa):
//            prior(startState),
//            likelihood(prior.nextSample(), pMakeObservation, pObserveIfPresent)
//    {
//        init();
//    }

    void init() {
        ConstrainedFactorisedDistribution<TRAJECTORY>::operator =(prior * likelihood);
        if(basis.basisVectors.size() == 0) {
            std::vector<double> entropiesByVar = prior.approximateEntropiesByVarIndex();
//            basis = Basis(entropiesByVar, *this, 6.0); // Entropy minimisation
            basis = Basis(*this); // norm-minimisation
            std::cout << "Basis entropy = " << basis.calculateEntropy(entropiesByVar) << std::endl;
        }
    }

    friend std::ostream &operator <<(std::ostream &out, const ABMPosterior<TRAJECTORY, STARTSTATE> &posterior) {
        out << posterior.prior;
        out << "Likelihood: " << posterior.likelihood;
        out << "Posteroir: " << static_cast<ConstrainedFactorisedDistribution<TRAJECTORY>>(posterior);
        out << "Basis (nVectors x domainSize): " << posterior.basis.basisVectors.size() << " x " << posterior.basis.origin.size() << std::endl;
        return out;
    }

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void load(Archive &ar, const unsigned int version) {
        ar >> prior >> likelihood >> basis;
        init();
    }

    template <typename Archive>
    void save(Archive &ar, const unsigned int version) const {
        ar << prior << likelihood << basis;
    }

//    template <typename Archive>
//    void serialize(Archive &ar, const unsigned int version) {
//        ar & pObserveIfPresent & pPredator & pPrey & kappa & realTrajectory & observations & basisObj;
//    }


    BOOST_SERIALIZATION_SPLIT_MEMBER();

};


template<class TRAJECTORY, class STARTSTATE>
ABMPosterior<TRAJECTORY,STARTSTATE> makeABMPosterior(STARTSTATE startState, double pMakeObservation, double pObserveIfPresent, double kappa) {
    ABMPrior<TRAJECTORY,STARTSTATE> prior(startState);
    ABMLikelihood<TRAJECTORY> likelihood(prior.nextSample(), pMakeObservation, pObserveIfPresent, kappa);
    return ABMPosterior<TRAJECTORY,STARTSTATE>(std::move(prior), std::move(likelihood));
}


#endif //ABMCMC_ABMPOSTERIOR_H
