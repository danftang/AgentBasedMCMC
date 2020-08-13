package hermiteForm

import kotlin.random.Random

class MCMCTreeNode(val value: IntArray, val breakpoint: Int, val parent: MCMCTreeNode?, val children: ArrayList<MCMCTreeNode> = ArrayList()) {

    val depth: Int
        get() = if(parent == null) 0 else parent.depth + 1


    fun splitAt(newBreakpoint: Int, dynamics: HermiteDecomposition): MCMCTreeNode {
        val reducedNullBasis = SparseColIntMatrix(dynamics.nullBasis.subList(0, newBreakpoint), dynamics.nullBasis.nRows)
        val fixedPartValue = value.mapIndexed { i, x ->
            if(i > newBreakpoint) x else if(i == newBreakpoint) 1-x else 0
        }
//        println("Trying with fixed basis $fixedPartValue")
        val baseActVector = dynamics.basisVectorToActs(fixedPartValue)
//        println("baseActVector is $baseActVector")
//        val constraintValues = dynamics.actToNonBasisVector(baseActVector)
        val newSolutionLowerPart = reducedNullBasis.IPsolve(-baseActVector, DoubleArray(reducedNullBasis.nCols) {0.0}.asList())
        val newSolution = IntArray(value.size) { i ->
            if(i < newSolutionLowerPart.size) newSolutionLowerPart[i] else fixedPartValue[i]
        }
        val child = MCMCTreeNode(newSolution, newBreakpoint, this)
        children.add(child)
        return child
    }

    fun generateAllChildren(decomposition: HermiteDecomposition): List<MCMCTreeNode> {
        for (n in 1 until breakpoint) {
            try {
                val newBreakpoint = n
                if (newBreakpoint != null) {
//                    println("Trying Break at $newBreakpoint")
                    val childSample = splitAt(newBreakpoint, decomposition)
                    val solution = decomposition.basisVectorToActs(childSample.value)
                    println("Sample: $solution")
//                    plotFeynmannDiagram(solution, abmMatrix, N)
                }
            } catch (err: java.lang.RuntimeException) {
//                    println("Got runtime exception")
            }
        }
        return children
    }

    fun chooseNextBreakpoint(): Int? {
        val newBreakpoint = breakpoint - Random.nextInt(0, breakpoint-1)
//        var newBreakpoint: Int
//        if(value.asList().subList(0, breakpoint).sum() == 0) return null
//        do {
//            newBreakpoint = breakpoint - Random.nextInt(0, breakpoint)
//        } while(value[newBreakpoint] == 0)
        return newBreakpoint
    }

}