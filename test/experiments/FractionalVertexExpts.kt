package experiments

import ABMCMC
import PredatorPreyABM
import org.junit.Test

class FractionalVertexExpts {

    @Test
    fun simpleFractionMCMC() {
        val predatorInitialDensity = 0.02
        val preyInitialDensity = 0.04
        val nTimesteps = 8
        PredatorPreyABM.gridSize = 16
        val (observations, realTrajectory) = PredatorPreyExpts.generateObservations(
            PredatorPreyABM.randomFermionicState(predatorInitialDensity, preyInitialDensity),
            nTimesteps,
            0.02
        )
        val prior = PredatorPreyABM.Prior(predatorInitialDensity, preyInitialDensity)

//        checkTrajectorySatisfiesObervations(realTrajectory, observations)
////        checkTrajectorySatisfiesObervationConstraints(realTrajectory, observations)
//        println("Checking real trajectory against observation constraints")
//        checkTrajectorySatisfiesConstraints(realTrajectory, observations.flatMap { it.eventConstraints() })
//        println("Checking real trajectory is fermionic")
//        checkTrajectorySatisfiesConstraints(realTrajectory, fermionicConstraints( nTimesteps, PredatorPreyABM))
//        println("Checking real trajectory is continuous")
//        checkTrajectorySatisfiesConstraints(realTrajectory, continuityConstraints( nTimesteps, PredatorPreyABM))

        val mcmc = ABMCMC(PredatorPreyABM, nTimesteps, observations + prior, realTrajectory.eventVector)
        println("Initial state is ${mcmc.simplex.X()}")
        println("Starting sampling")

        for(s in 1..5000) {
            val sample = mcmc.fractionSample()
        }
    }
}