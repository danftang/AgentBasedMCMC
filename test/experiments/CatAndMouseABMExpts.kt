package experiments

import ABMCMC
import CatAndMouseABM
import Trajectory
import org.junit.Test

class CatAndMouseABMExpts {

    @Test
    fun fermionic2Timestep() {
        val observations = listOf(CatAndMouseABM.CMObservation(
            CatAndMouseABM.CatMouseAgent(CatAndMouseABM.AgentType.CAT,CatAndMouseABM.AgentPosition.LEFT),
            1,
            true
        ))
        val mcmc = ABMCMC(CatAndMouseABM, 2, observations)
        assert(mcmc.simplex.isPrimalFeasible())

        val initialTrajectory = Trajectory(CatAndMouseABM, mcmc.simplex.X())
        println("Initial Trajectory")
        println(initialTrajectory)
        CatAndMouseABM.plot(initialTrajectory)

        for(sample in 1..6) {
            val sample = mcmc.nextSample()
            CatAndMouseABM.plot(sample)
            println(sample)
        }
    }
}