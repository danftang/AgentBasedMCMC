package experiments

import ABMCMC
import MinimisationMCMC
import PredatorPreyABM
import Simplex
import SimplifiedSimplex
import org.junit.Test
import kotlin.math.exp
import kotlin.math.roundToInt
import kotlin.random.Random

class uniformDegeneracyProbExpts {

    // Do a random walk around vertices and watch saprsity
    @Test
    fun randomWalkSparsity() {
        val predatorInitialDensity = 0.02
        val preyInitialDensity = 0.04
        val nTimesteps = 4
        PredatorPreyABM.gridSize = 16
        val (observations, realTrajectory) = PredatorPreyExpts.generateObservations(
            PredatorPreyABM.randomFermionicState(predatorInitialDensity, preyInitialDensity),
            nTimesteps,
            0.02
        )
        val prior = PredatorPreyABM.Prior(predatorInitialDensity, preyInitialDensity)

        val mcmc = ABMCMC(PredatorPreyABM, nTimesteps, observations + prior, realTrajectory.eventVector)
        println("Initial state is ${mcmc.simplex.X()}")
        println("Starting sampling")

        for(s in 1..5000) {
            if(s.rem(10) == 0) {
//                println("$s sparsity = ${mcmc.simplex.M.sparsity()} size = ")
                println("$s")
            }
            mcmc.simplex.randomPivot()
        }
    }


    @Test
    fun simplifiedRandomWalkSparsity() {
        val predatorInitialDensity = 0.02
        val preyInitialDensity = 0.04
        val nTimesteps = 4
        PredatorPreyABM.gridSize = 16
        val (observations, realTrajectory) = PredatorPreyExpts.generateObservations(
            PredatorPreyABM.randomFermionicState(predatorInitialDensity, preyInitialDensity),
            nTimesteps,
            0.02
        )
        val prior = PredatorPreyABM.Prior(predatorInitialDensity, preyInitialDensity)
        val constraints = ABMCMC.constraints(PredatorPreyABM, nTimesteps, observations)

        val mcmc = SimplifiedSimplex(constraints, emptyMap(), realTrajectory.eventVector.nonZeroEntries)
        println("Initial state is ${mcmc.X()}")
        println("Starting sampling")

        for(s in 1..5000) {
            if(s.rem(10) == 0) {
                println("$s size = ${mcmc.columns.sumBy { it.size } } / ${mcmc.columns.size * mcmc.rows.size}")
//                println("$s")
            }
            mcmc.randomPivot()
        }
    }


    @Test
    fun minimisationMCMCTest() {
        val predatorInitialDensity = 0.02
        val preyInitialDensity = 0.04
        val nTimesteps = 4
        PredatorPreyABM.gridSize = 16
        val (observations, realTrajectory) = PredatorPreyExpts.generateObservations(
            PredatorPreyABM.randomFermionicState(predatorInitialDensity, preyInitialDensity),
            nTimesteps,
            0.02
        )
        val prior = PredatorPreyABM.Prior(predatorInitialDensity, preyInitialDensity)
        val constraints = ABMCMC.constraints(PredatorPreyABM, nTimesteps, observations)

//        val mcmc = MinimisationMCMC(constraints, realTrajectory.eventVector.nonZeroEntries)
        val mcmc = MinimisationMCMC(constraints)
        println("Problem size = ${mcmc.nConstraints} x ${mcmc.nVariables}")

        for(s in 1..10) {
            val sample = mcmc.nextSample()
            println("objective is ${mcmc.solver.objective().value()}")
//            println(sample)
        }
//        println("primes = ${MinimisationMCMC.logPrimes.map { exp(it).roundToInt() }}")
    }
}

fun<T> Simplex<T>.randomPivot() where T : Comparable<T>, T: Number {
    var pivotCol: Int
    var pivotRows: List<Int>
    do {
        pivotCol = Random.nextInt(M.nCols-1)
        pivotRows = pivotableRows(pivotCol, false)
    } while(pivotRows.isEmpty())
    pivot(pivotRows.random(),pivotCol)
}

fun SimplifiedSimplex.randomPivot() {
    var pivotCol: Int
    var pivotRows: List<Int>
    do {
        pivotCol = Random.nextInt(nVariables)
        pivotRows = pivotableRows(pivotCol, false)
    } while(pivotRows.isEmpty())
    pivot(pivotRows.random(),pivotCol)
}