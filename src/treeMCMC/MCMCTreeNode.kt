package treeMCMC

import lib.SparseColIntMatrix
import org.apache.commons.math3.distribution.EnumeratedDistribution
import org.apache.commons.math3.util.Pair
import java.lang.RuntimeException
import kotlin.math.exp
import kotlin.math.min
import kotlin.math.pow
import kotlin.random.Random

// breakpoint is the basis that is definitely different from parent, and takes on value breakpointValue
// all bases after breakpoint are identical to parent
// all bases before breakpoint are recalculated with IPSolve
class MCMCTreeNode(val tree: MCMCTree, val parent: MCMCTreeNode?, val breakpoint: Int, val breakpointValue: Int, var pTransitionFromParent: Double) {

    val basisSize: Int
        get() = tree.hermite.nullBasis.nCols

    // If p is parent prob, c is child prob, t_pc is transition prob from parent to child
    // and t_cp is transition probability from child to parent
    // then acceptance probability of transition from parent to child is
    // Pa = max(1, c/p . t_cp/t_pc)
    // so if t_cp = t_pc . p/c then acceptance is 1
    // To ensure we have some probability left over for children we cap this at P_cap
    // so t_cp = min(P_cap, t_pc . p/c)
    val pTransitionToParent: Double
        get() = min(Pcap, pTransitionFromParent * exp(parent?.targetLnProbability?:Double.NEGATIVE_INFINITY - targetLnProbability))

    val pAcceptTransitionFromParent: Double
        get() {
            try {
                return if(parent != null) {
                    Pcap / (pTransitionFromParent * exp(parent.targetLnProbability - targetLnProbability))
                } else {
                    0.0
                }
            } catch(err: RuntimeException) {
                println("Infeasible child")
                return 0.0
            }
        }

    val actValue: SparseColIntMatrix.SparseIntColumn
        get() = tree.hermite.basisVectorToActs(basisValue)

    // Lazy values
    val targetLnProbability: Double by lazy {
        basisValue.asSequence().withIndex().sumByDouble {it.value * tree.basisLnProbabilities[it.index] }
    }

    val basisValue: IntArray by lazy {
        val reducedNullBasis = SparseColIntMatrix(tree.hermite.nullBasis.subList(0, breakpoint), tree.hermite.nullBasis.nRows)
        val newValue = IntArray(basisSize) {0}
        if(parent != null) {
            for(i in breakpoint+1 until basisSize) {
                newValue[i] = parent.basisValue[i]
            }
            newValue[breakpoint] = breakpointValue
        }
        val baseActVector = tree.hermite.basisVectorToActs(newValue)
        val newSolutionLowerPart = reducedNullBasis.IPsolve(-baseActVector, DoubleArray(reducedNullBasis.nCols) {0.0}.asList())
        for(i in 0 until breakpoint) {
            newValue[i] = newSolutionLowerPart[i]
        }
        newValue
    }

    val children: EnumeratedDistribution<MCMCTreeNode> by lazy {
        val acts = SparseColIntMatrix.SparseIntColumn()

        // calculate all upper bounds
        val upperBounds = basisUpperBounds()

        val childNodes = ArrayList<MCMCTreeNode>()
        val currentChildTransitionProb = 1.0 - pTransitionToParent
        for(col in 0 until breakpoint) {
            val upperBound = upperBounds[col]
            if(upperBound > 0) {
                for(breakpointVal in 0..upperBound) {
                    if(breakpointVal != basisValue[col]) {
                        childNodes.add(MCMCTreeNode(tree, this, col, breakpointVal, tree.actProbabilities[col].pow(breakpointVal)))
                    }
                }
            }
        }
        val childTotalProb = childNodes.sumByDouble { it.pTransitionFromParent }
        childNodes.forEach { child ->
            child.pTransitionFromParent  *= currentChildTransitionProb/childTotalProb
        }
        EnumeratedDistribution(childNodes.map { Pair(it, it.pTransitionFromParent) })
    }



    fun basisUpperBounds(): IntArray {
        val acts = SparseColIntMatrix.SparseIntColumn()
        val upperBounds = IntArray(breakpoint) {0} // Inclusive upper bound for each basis value up to breakpoint
        for(col in basisSize-1 downTo 0) {
            val thisBasis = tree.hermite.nullBasis[col]
            if(col < breakpoint) {
                thisBasis.entries.filter { it.value < 0 }.forEach { antiAgent ->
                    if (tree.lastPositiveBasis[antiAgent.key] ?: -1 < col)
                        upperBounds[col] = (tree.hermite.baseState[antiAgent.key] - acts[antiAgent.key]) / antiAgent.value
                }
            }
            acts.weightedPlusAssign(thisBasis, basisValue[col])
        }
        return upperBounds
    }


    fun getNextSample(): MCMCTreeNode {
        if (parent != null && Random.nextDouble() < pTransitionToParent) {
            return parent // always accept transition to parent
        } else {
            val proposal = children.sample()
            if(Random.nextDouble() < proposal.pAcceptTransitionFromParent) return proposal
        }
        return this
    }


    companion object {
        const val Pcap = 0.5
    }


//    fun fixedBasisValue(index: Int): Int {
//        if(parent == null) return 0
//        return if(index > breakpoint) parent.value[index] else if(index == breakpoint) breakpointValue else 0
//    }


//    fun calculatePreBreakpointSolution(dynamics: HermiteDecomposition, basisLnProbabilities: DoubleArray) {
//        val reducedNullBasis = SparseColIntMatrix(dynamics.nullBasis.subList(0, breakpoint), dynamics.nullBasis.nRows)
//        val fixedPartValue = value.mapIndexed { i, x ->
//            if(i >= breakpoint) x else 0
//        }
//        val baseActVector = dynamics.basisVectorToActs(fixedPartValue)
//        val newSolutionLowerPart = reducedNullBasis.IPsolve(-baseActVector, DoubleArray(reducedNullBasis.nCols) {0.0}.asList())
//        for(i in 0 until breakpoint) {
//            value[i] = newSolutionLowerPart[i]
//        }
//    }


//    private fun targetLnProbability(basisLnProbabilities: DoubleArray): Double {
//        return value.asSequence().withIndex().sumByDouble {it.value * basisLnProbabilities[it.index] }
//    }





//    fun splitAt(newBreakpoint: Int, dynamics: HermiteDecomposition): MCMCTreeNode {
//        val reducedNullBasis = SparseColIntMatrix(dynamics.nullBasis.subList(0, newBreakpoint), dynamics.nullBasis.nRows)
//        val fixedPartValue = value.mapIndexed { i, x ->
//            if(i > newBreakpoint) x else if(i == newBreakpoint) 1-x else 0
//        }
////        println("Trying with fixed basis $fixedPartValue")
//        val baseActVector = dynamics.basisVectorToActs(fixedPartValue)
////        println("baseActVector is $baseActVector")
////        val constraintValues = dynamics.actToNonBasisVector(baseActVector)
//        val newSolutionLowerPart = reducedNullBasis.IPsolve(-baseActVector, DoubleArray(reducedNullBasis.nCols) {0.0}.asList())
//        val newSolution = IntArray(value.size) { i ->
//            if(i < newSolutionLowerPart.size) newSolutionLowerPart[i] else fixedPartValue[i]
//        }
//        val child = MCMCTreeNode(newSolution, newBreakpoint, this)
//        children.add(child)
//        return child
//    }

//    fun generateAllChildren(decomposition: HermiteDecomposition) {
//        for (n in 1 until breakpoint) {
//            try {
//                val newBreakpoint = n
//                if (newBreakpoint != null) {
////                    println("Trying Break at $newBreakpoint")
//                    val childSample = splitAt(newBreakpoint, decomposition)
//                    val solution = decomposition.basisVectorToActs(childSample.value)
//                    println("Sample: $solution")
////                    plotFeynmannDiagram(solution, abmMatrix, N)
//                }
//            } catch (err: java.lang.RuntimeException) {
////                    println("Got runtime exception")
//            }
//        }
//    }

//    fun chooseNextBreakpoint(): Int? {
//        val newBreakpoint = breakpoint - Random.nextInt(0, breakpoint-1)
////        var newBreakpoint: Int
////        if(value.asList().subList(0, breakpoint).sum() == 0) return null
////        do {
////            newBreakpoint = breakpoint - Random.nextInt(0, breakpoint)
////        } while(value[newBreakpoint] == 0)
//        return newBreakpoint
//    }


}