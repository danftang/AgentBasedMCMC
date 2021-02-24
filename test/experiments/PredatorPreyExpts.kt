package experiments

import ABMCMC
import Constraint
import PredatorPreyABM
import lib.collections.Multiset
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.random.Random

class PredatorPreyExpts {

    @Test
    fun fermionicPredPrey() {
        val nTimesteps = 4
        PredatorPreyABM.gridSize = 4
//        val observations = generateObservations(
//            PredatorPreyABM.randomState(0.2, 0.3),
//            nTimesteps,
//            0.02
//        )
//        val observationConstraints = observations.flatMap { it.constraints() }

        val observationConstraints = listOf<Constraint<Fraction>>()

        val mcmc = ABMCMC(PredatorPreyABM, nTimesteps, observationConstraints)
        for(sample in 1..10) {
            val sample = mcmc.nextSample()
            //println(mcmc.nextSample())
        }

    }

    fun generateObservations(
        startState: Multiset<PredatorPreyABM.PredPreyAgent>,
        nTimesteps: Int,
        pMakeObservation: Double): List<PredatorPreyABM.CompletedObservation> {
        val trajectory = PredatorPreyABM.runABM(startState, nTimesteps)
        val observations = ArrayList<PredatorPreyABM.CompletedObservation>(
            (nTimesteps*PredatorPreyABM.agentDomain.size*pMakeObservation).toInt()
        )
        for(t in 0 until nTimesteps) {
            for(x in 0 until PredatorPreyABM.gridSize) {
                for(y in 0 until PredatorPreyABM.gridSize) {
                    if(Random.nextDouble() < pMakeObservation) {
                        observations.add(
                            PredatorPreyABM.CompletedObservation(
                                PredatorPreyABM.PredPreyAgent(x,y,PredatorPreyABM.AgentType.PREDATOR), t, trajectory)
                        )
                        observations.add(
                            PredatorPreyABM.CompletedObservation(
                                PredatorPreyABM.PredPreyAgent(x,y,PredatorPreyABM.AgentType.PREY), t, trajectory)
                        )
                    }
                }
            }
        }
        return observations
    }


}