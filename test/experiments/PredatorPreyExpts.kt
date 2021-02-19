package experiments

import ABMCMC
import Constraint
import PredatorPreyABM
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test

class PredatorPreyExpts {

    @Test
    fun fermionicPredPrey() {
        val timesteps = 8
        val observations = listOf<Constraint<Fraction>>()
        val mcmc = ABMCMC(PredatorPreyABM, timesteps, observations)

        for(sample in 1..10) {
            val sample = mcmc.nextSample()
            //println(mcmc.nextSample())
        }

    }
}