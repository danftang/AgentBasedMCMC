package experiments

import ABMCMC
import ABMCMC.Companion.validTrajectoryConstraints
import MutableConstraint
import PredatorPreyABM
import Trajectory
import isSatisfiedBy
import lib.abstractAlgebra.FractionOperators
import lib.collections.Multiset
import lib.sparseVector.emptySparseVector
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.random.Random

class PredatorPreyExpts {

    @Test
    fun fermionicPredPrey() {
        val predatorInitialDensity = 0.02
        val preyInitialDensity = 0.04
        val nTimesteps = 8
        PredatorPreyABM.gridSize = 16
        val (observations, realTrajectory) = generateObservations(
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

        val expectation = mcmc.expectation(10000, emptySparseVector(FractionOperators)) { sample, sum ->
            sample + sum
        }
        val expectationTrajectory = Trajectory(PredatorPreyABM, expectation)
//        for(n in 1..10000) {
//            val sample = mcmc.nextSample()
//            if(!mcmc.simplex.isInteger()) println("Fractional solution. Log fraction penalty = ${mcmc.simplex.logFractionPenalty(mcmc.simplex.X())}")
//            histogram.add(sample)
//            //println(mcmc.nextSample())
//        }

        val plotTime = nTimesteps
        println("Histogram at time $plotTime ${expectationTrajectory.stateAt(plotTime)}")
        with(PredatorPreyABM) {
            plotHeatMap(expectationTrajectory.stateAt(plotTime))
                .replotPoints(realTrajectory.stateAt(plotTime))
                .renderAndClose()
        }
    }


    fun testDoubleSimplex() {

    }


    @Test
    fun testORSolve() {
        val nTimesteps = 8
        PredatorPreyABM.gridSize = 20
        val (observations, _) = generateObservations(
            PredatorPreyABM.randomFermionicState(0.2, 0.3),
            nTimesteps,
            0.02
        )
        val constraints = validTrajectoryConstraints(PredatorPreyABM, nTimesteps,false) + observations.flatMap { it.eventConstraints() }
        val objective =
            //observations.flatMap { observation -> observation.eventConstraints().flatMap { it.coefficients.keys } }.associateWith { Fraction.ONE }
            //(0 until constraints.numVars()).associateWith { Fraction.ONE }
            emptyMap<Int,Fraction>()

//        val simplex = Simplex(constraints, objective.asVector(FractionOperators))
//        simplex.pivotToInitialIntegerSolution()
//        println("Initial solution is ${simplex.X()}")

        val solution = ORTools.BooleanSolve(constraints, objective)
        println("solution is ${solution.asList()}")
    }

    companion object {

        fun generateObservations(
            startState: Multiset<PredatorPreyABM.PredPreyAgent>,
            nTimesteps: Int,
            pMakeObservation: Double
        ): Pair<List<PredatorPreyABM.PPObservation>, Trajectory<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts>> {
            val trajectory = PredatorPreyABM.fermionicRunABM(startState, nTimesteps)
            val observations = ArrayList<PredatorPreyABM.PPObservation>(
                (nTimesteps * PredatorPreyABM.agentDomain.size * pMakeObservation).toInt()
            )
            for (t in 0 until nTimesteps) {
                for (x in 0 until PredatorPreyABM.gridSize) {
                    for (y in 0 until PredatorPreyABM.gridSize) {
                        if (Random.nextDouble() < pMakeObservation) {
                            observations.add(
                                PredatorPreyABM.PPObservation(
                                    t,
                                    PredatorPreyABM.PredPreyAgent(x, y, PredatorPreyABM.AgentType.PREDATOR),
                                    trajectory
                                )
                            )
                            observations.add(
                                PredatorPreyABM.PPObservation(
                                    t, PredatorPreyABM.PredPreyAgent(x, y, PredatorPreyABM.AgentType.PREY), trajectory
                                )
                            )
                        }
                    }
                }
            }
            checkTrajectorySatisfiesObervations(trajectory, observations)
            return Pair(observations, trajectory)
        }

        fun checkTrajectorySatisfiesObervations(
            trajectory: Trajectory<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts>,
            observations: List<PredatorPreyABM.PPObservation>
        ) {
            assert(observations.all { it.logLikelihood(trajectory) != Double.NEGATIVE_INFINITY })
        }


        fun checkTrajectorySatisfiesObervationConstraints(
            trajectory: Trajectory<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts>,
            observations: List<PredatorPreyABM.PPObservation>
        ) {
            assert(observations.all { it.eventConstraints().all { it.isSatisfiedBy(trajectory.eventVector) } })
        }


        fun checkTrajectorySatisfiesConstraints(
            trajectory: Trajectory<PredatorPreyABM.PredPreyAgent, PredatorPreyABM.Acts>,
            constraints: List<MutableConstraint<Fraction>>
        ) {
            assert(constraints.all { it.isSatisfiedBy(trajectory.eventVector) })
        }

    }
}