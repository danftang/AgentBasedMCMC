package experiments

import ABMCMC
import IntegerSimplexMCMC
import PredatorPreyABM
import lib.sparseVector.SparseVector
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.math.exp
import kotlin.math.min
import kotlin.random.Random

class FractionalVertexExpts {

    @Test
    fun simpleFractionMCMC() {
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
            if(s.rem(100) == 0) println("Sample $s")
            val sample = mcmc.simplex.fractionSampleV2()
        }
    }
}


fun IntegerSimplexMCMC<Fraction>.fractionSample(): SparseVector<Fraction> {
    val proposal = proposePivot()

    nSamples++
    if(state.logFractionalPenalty == 0.0) {
        forwardPivot(proposal)
        if(state.logFractionalPenalty != 0.0) {
            println("$nSamples penalty = ${state.logFractionalPenalty}")
        }
        return X(false)
    }
    nFractionalSamples++
    forwardPivot(proposal)
    val acceptance = exp(min(0.0, state.logFractionalPenalty - revertState.cache.logFractionalPenalty ))
    if(Random.nextDouble() <= acceptance) {
        println("$nSamples penalty = ${state.logFractionalPenalty}")
        if(state.logFractionalPenalty == 0.0) {
            println("Fractional runlength was $nFractionalSamples samples")
            nFractionalSamples = 0
        }
    } else {
        println("rejecting")
        revertLastPivot()
    }
    return X(false)
}

fun IntegerSimplexMCMC<Fraction>.fractionSampleV2(): SparseVector<Fraction> {
    val proposal = proposePivot()

    nSamples++
    if(state.logFractionalPenalty == 0.0) {
        nFractionalSamples = 0
        forwardPivot(proposal)
        if(state.logFractionalPenalty != 0.0) {
            println("$nSamples penalty = ${state.logFractionalPenalty}")
        }
        return X(false)
    }
    nFractionalSamples++

    val proposedLogFractionalPenalty = logFractionalPenaltyAfterPivot(proposal.col)
    val logAcceptance = min(0.0,proposedLogFractionalPenalty - state.logFractionalPenalty)
    if(Random.nextDouble() < exp(logAcceptance)) {
        println("$nFractionalSamples Accepted. Penalty = $proposedLogFractionalPenalty")
        forwardPivot(proposal, proposedLogFractionalPenalty)
    } else {
        println("$nFractionalSamples Rejected. Penalty = $proposedLogFractionalPenalty")
    }
    return X(false)
}

