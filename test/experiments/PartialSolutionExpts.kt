package experiments

import ABMMatrices.twoDabmMatrix
import ConvexPolyhedron
import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.HashRowColIntMatrix
import org.junit.Test
import kotlin.random.Random

// Propose constraints on the solution by specifying
// a subset of act values and using MIP solver to
// find a completion or prove infeasibility
class PartialSolutionExpts {

    init {
        System.loadLibrary("jniortools")
    }


    // Find out how often we get infeasible solutions
    // if we constrain on random subsets of acts
    @Test
    fun infieasibilityProportion() {
        val gridSize = 32
        val timesteps = 6
        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
        println("abm sparsity is ${abmMatrix.sparsityRatio()}")
        println("abm geometry is ${abmMatrix.nRows} x ${abmMatrix.nCols}")
//        println(abmMatrix.toSparsityString())
        val observations = HashIntVector()
        for(agent in 1..20) {
            val xPos = Random.nextInt(gridSize)
            val yPos = Random.nextInt(gridSize)
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }

        // zero removal
        val polyhedron = ConvexPolyhedron(HashRowColIntMatrix(abmMatrix), HashIntVector(observations))
        println("Original polyhedron size is ${polyhedron.M.nRows} ${polyhedron.M.nCols}")
        polyhedron.removeZeros()
        println("Reduced polyhedron size is ${polyhedron.M.nRows} ${polyhedron.M.nCols}")

        for(nConstraints in 0 until 100) {
            println("Starting solve with $nConstraints constraints")
            val solution = polyhedron.findValidPoint()
            //println(SparseColIntMatrix.SparseIntColumn(solution))
            println("Error vector (should be empty) = ${abmMatrix * solution - observations}")
            val newConstraint = HashIntVector(Random.nextInt(polyhedron.M.nCols) to 1)
            polyhedron.constrainSolution(newConstraint)
        }
    }
}