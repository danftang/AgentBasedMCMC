package experiments

import ABMCMC
import CatAndMouseABM
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

        for(sample in 1..3) {
            val sample = mcmc.nextSample()
            CatAndMouseABM.plot(sample)
            println(sample)
        }
    }
}