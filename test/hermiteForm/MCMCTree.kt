package hermiteForm

import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.util.Pair
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.pow
import kotlin.random.Random

class MCMCTree {
    var currentSample: MCMCTreeNode
    val abmMatrix: SparseColIntMatrix
    val hermite: HermiteDecomposition
    val actProbabilities: DoubleArray
    val basisLnProbabilities: DoubleArray
    val lastPositiveBasis: Map<Int,Int>  // map from an act to the left-most basis with a positive coefficient in that act

    constructor(abmMatrix: SparseColIntMatrix, observations: SparseColIntMatrix.SparseIntColumn, actProbabilities: DoubleArray) {
        this.actProbabilities = actProbabilities
        this.abmMatrix = abmMatrix
        hermite = HermiteDecomposition(abmMatrix, observations)
        lastPositiveBasis = HashMap()
        for(col in hermite.nullBasis.nCols-1 downTo 0) {
            for(act in hermite.nullBasis[col].entries) {
                if(act.value > 0) lastPositiveBasis[act.key] = col
            }
        }
        basisLnProbabilities = hermite.nullBasis.map { basis ->
            basis.entries.sumByDouble { it.value * ln(actProbabilities[it.key]) }
        }.toDoubleArray()
        currentSample = MCMCTreeNode(this, null, hermite.nullBasis.nCols, -1, 0.0)

    }

    fun nextSample() {
        currentSample = currentSample.getNextSample()
    }

//    fun getNextSample(): IntArray {
//        var pAcceptNewSample: Double
//        var newSample: MCMCTreeNode
//        do {
//            val pChildParentTransition = currentSample.childParentTransitionProb()
//            if (Random.nextDouble() < pChildParentTransition) { // transition to parent
//                newSample = currentSample.parent!!
//                pAcceptNewSample = currentSample.parentChildTransitionProb / pChildParentTransition
//            } else {
//                if (currentSample.children == null) { // calculate children
//                    currentSample.children = allChildProposals()
//                }
//                newSample = currentSample.children!!.sample()
//                currentSample.calculatePreBreakpointSolution(hermite)
//                pAcceptNewSample = newSample.childParentTransitionProb()/newSample.parentChildTransitionProb
//            }
//            pAcceptNewSample = exp(newSample.targetLnProbability(basisLnProbabilities) - currentSample.targetLnProbability(basisLnProbabilities))
//        } while(Random.nextDouble() < pAcceptNewSample)
//        return currentSample.value
//    }


    // gets all children that aren't immediately discounted
    // by the null vectors (i.e. any basis with an anti-agent in a position where there are no bases
    // to the left with agents in that position - left of the left-most positive coefficient -
    // must not drive the agent count negative)
//    fun allChildProposals(): EnumeratedDistribution<MCMCTreeNode> {
//        val acts = SparseColIntMatrix.SparseIntColumn()
//        val upperBounds = HashMap<Int,Int>() // basis to inclusive upper bound
//        for(col in currentSample.value.size-1 downTo 0) {
//            val thisBasis = hermite.nullBasis[col]
//            thisBasis.entries.filter { it.value < 0 }.forEach { antiAgent ->
//                if(lastPositiveBasis[antiAgent.key]?:-1 < col) upperBounds[col] = (hermite.baseState[antiAgent.key] - acts[antiAgent.key])/antiAgent.value
//            }
//            acts.weightedPlusAssign(thisBasis, currentSample.value[col])
//        }
//        val childPMF = ArrayList<org.apache.commons.math3.util.Pair<MCMCTreeNode,Double>>()
//        val currentChildTransitionProb = 1.0-currentSample.childParentTransitionProb(basisLnProbabilities)
//        for(col in 0 until currentSample.value.size) {
//            val upperBound = upperBounds[col]?:0
//            if(upperBound > 0) {
//                for(v in 0..upperBound) {
//                    if(v != currentSample.value[col]) {
//                        val newSample = currentSample.value.copyOf()
//                        newSample[col] = v
//                        childPMF.add(
//                            Pair(
//                                MCMCTreeNode(
//                                    newSample,
//                                    col,
//                                    currentSample,
//                                    actProbabilities[col]
//                                ), actProbabilities[col].pow(v)
//                            )
//                        )
//                    }
//                }
//            }
//        }
//        val childTotalProb = childPMF.sumByDouble { it.second }
//        childPMF.forEach { pair ->
//            pair.first.parentChildTransitionProb  *= currentChildTransitionProb/childTotalProb
//        }
//        return EnumeratedDistribution(childPMF)
//    }
}