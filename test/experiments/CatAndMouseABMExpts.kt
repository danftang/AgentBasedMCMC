package experiments

import ABMCMC
import CatAndMouseABM
import org.junit.Test

class CatAndMouseABMExpts {

    @Test
    fun fermionic2Timestep() {
        val mcmc = ABMCMC(CatAndMouseABM, 2, emptyList())

        for(sample in 1..10) {
            println(mcmc.nextSample())
        }
    }
}