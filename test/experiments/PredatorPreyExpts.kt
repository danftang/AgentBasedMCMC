package experiments

import ABMCMC
import ABMCMC.Companion.validTrajectoryConstraints
import Constraint
import PredatorPreyABM
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test

class PredatorPreyExpts {

    @Test
    fun fermionicPredPrey() {
        val timesteps = 8
        PredatorPreyABM.gridSize = 8
        val observations = listOf<Constraint<Fraction>>(
            Constraint(hashMapOf(16 to Fraction.ONE), "==", Fraction.ONE),
            Constraint(hashMapOf(60 to Fraction.ONE), "==", Fraction.ONE)
        )

//        val startState = ORTools.GlopSolve(validTrajectoryConstraints(PredatorPreyABM, timesteps) + observations)
//        println(startState.size)
//        for(i in 0 until startState.size) {
//            if(startState[i] != 0.0) println("$i -> ${startState[i]}")
//        }

        val mcmc = ABMCMC(PredatorPreyABM, timesteps, observations)
        for(sample in 1..10) {
            val sample = mcmc.nextSample()
            //println(mcmc.nextSample())
        }

    }
}