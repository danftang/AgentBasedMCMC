package experiments

import ABMCMC
import PredatorPreyABM
import Trajectory
import lib.collections.Multiset
import org.junit.Test
import kotlin.random.Random

class PredatorPreyExpts {

    @Test
    fun fermionicPredPrey() {
        val nTimesteps = 8
        PredatorPreyABM.gridSize = 32
        val observations = generateObservations(
            PredatorPreyABM.randomState(0.2, 0.3),
            nTimesteps,
            0.02
        )

        val mcmc = ABMCMC(PredatorPreyABM, nTimesteps, observations)
        println("Initial state is ${mcmc.simplex.X()}")
        println("Starting sampling")
        for(sample in 1..1000) {
            val sample = mcmc.nextSample()
            if(!mcmc.simplex.isInteger()) println("Fractional solution. Log fraction penalty = ${mcmc.simplex.logFractionPenalty(mcmc.simplex.X())}")
            //println(mcmc.nextSample())
        }

    }

    fun generateObservations(
        startState: Multiset<PredatorPreyABM.PredPreyAgent>,
        nTimesteps: Int,
        pMakeObservation: Double): List<PredatorPreyABM.PPObservation> {
        val trajectory = PredatorPreyABM.runABM(startState, nTimesteps)
        val observations = ArrayList<PredatorPreyABM.PPObservation>(
            (nTimesteps*PredatorPreyABM.agentDomain.size*pMakeObservation).toInt()
        )
        for(t in 0 until nTimesteps) {
            for(x in 0 until PredatorPreyABM.gridSize) {
                for(y in 0 until PredatorPreyABM.gridSize) {
                    if(Random.nextDouble() < pMakeObservation) {
                        observations.add(
                            PredatorPreyABM.PPObservation(
                                t, PredatorPreyABM.PredPreyAgent(x,y,PredatorPreyABM.AgentType.PREDATOR), trajectory
                            )
                        )
                        observations.add(
                            PredatorPreyABM.PPObservation(
                                t, PredatorPreyABM.PredPreyAgent(x,y,PredatorPreyABM.AgentType.PREY), trajectory
                            )
                        )
                    }
                }
            }
        }
        checkTrajectorySatisfiesObervations(trajectory, observations)
        return observations
    }

    fun checkTrajectorySatisfiesObervations(trajectory: Trajectory<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts>,
    observations: List<PredatorPreyABM.PPObservation>) {
        assert(observations.all { it.logLikelihood(trajectory) != Double.NEGATIVE_INFINITY })
    }


}