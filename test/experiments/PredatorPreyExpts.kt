package experiments

import ABMCMC
import Constraint
import PredatorPreyABM
import org.apache.commons.math3.fraction.Fraction

class PredatorPreyExpts {

    fun fermionicPredPrey() {
        val timesteps = 2
        val observations = listOf<Constraint<Fraction>>()
        val mcmc = ABMCMC(PredatorPreyABM, timesteps, observations)
    }
}