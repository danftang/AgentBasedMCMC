package experiments

import ABMMatrices.twoDabmMatrix
import Simplex
import lib.abstractAlgebra.DoubleOperators
import lib.collections.GridMap
import lib.sparseMatrix.GridMapMatrix
import org.junit.Test

// Experiments with pivoting the ABM matrix, a-la Simplex algorithm
class SimplexExpts {

    init {
        System.loadLibrary("jniortools")
    }

    @Test
    fun introToLinearProgrammingAndGameTheoryPage89() {
        val coeffs = GridMapMatrix(DoubleOperators,4,8)
        coeffs[0,0] = -2.0
        coeffs[0,2] = 6.0
        coeffs[0,3] = 2.0
        coeffs[0,5] = -3.0
        coeffs[0,6] = 1.0

        coeffs[1,0] = -4.0
        coeffs[1,1] = 1.0
        coeffs[1,2] = 7.0
        coeffs[1,3] = 1.0
        coeffs[1,5] = -1.0

        coeffs[2,2] = -5.0
        coeffs[2,3] = 3.0
        coeffs[2,4] = 1.0
        coeffs[2,5] = -1.0

        // constants
        coeffs[0,7] = 20.0
        coeffs[1,7] = 10.0
        coeffs[2,7] = 60.0

        // objective
        coeffs[3,2] = 13.0
        coeffs[3,3] = -6.0
        coeffs[3,5] = 2.0

        val initialPivots = intArrayOf(6,1,4)
        val simplex = Simplex(coeffs, initialPivots)
        println(simplex.M)
        simplex.minimise()
        println()
        println(simplex.X())
    }


//    @Test
//    fun oneDABM() {
//        val timesteps = 5
//        val N = 1
//        val abmMatrix = ABMMatrices.oneDNoInteractionMatrix(timesteps, N)
//        val observations = HashIntVector()
//        observations[0] = -1
//        observations[timesteps*N] = 1
//
////        println(abmMatrix)
//        val simplex = Simplex(abmMatrix, observations)
//        println(simplex)
//        simplex.pivotRootedSourceTree()
////        simplex.pivotAt(4,3)
////        println(simplex)
////        simplex.pivotAt(3,2)
////        println(simplex)
////        simplex.pivotAt(2,1)
////        println(simplex)
////        simplex.pivotAt(1,0)
//        println(simplex)
//        println(simplex.findAllPositivePivotCols())
//    }

//    @Test
//    fun twoDABM() {
//        val gridSize = 8
//        val timesteps = 4
//        val nAgents = 32
//        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
//        val observations = HashIntVector()
//        val startPos = Array(nAgents) {
//            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
//        }
//        for(agent in startPos.indices) {
//            val xPos = startPos[agent].first
//            val yPos = startPos[agent].second
//            observations[xPos*gridSize + yPos] = -1
//            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
//        }
//
//        val eventProbs = (0 until abmMatrix.nCols).associateWith { 1.0/abmMatrix.nCols }
//
//        val simplex = Simplex(abmMatrix, observations, eventProbs)
//
//        val Xinit = constructInitialSolution(simplex, startPos, gridSize, timesteps)
//
//        println("Starting Hermite transformation")
//        val hermiteTime = measureTimeMillis {
////            simplex.pivotRootedSourceTree(Xinit)
//            simplex.pivotSourcePolyTree(Xinit)
//        }
//        println("Done in $hermiteTime ms")
//
//        println("pivoted rows = ${simplex.pivotPoints.size}")
//        println("Total M size = ${simplex.M.nRows} x ${simplex.M.nCols}")
//        println("Base solution = ${simplex.X}")
////        println("Observation error = ${simplex.observationNorm()}")
//
//        assert(simplex.M*simplex.X == simplex.B)
//        assert(simplex.X.isPositive())
//
//        println("Row sizes = ${simplex.M.rows.map { it.sparseSize }}")
//        println("Col sizes = ${simplex.M.columns.map { it.sparseSize }}")
////        println("Smallest root reduction column ${
////        simplex.M.columns.filter {
////            it[simplex.M.nRows-1] != 0
////        }.map { it.sparseSize }.run { min() }
////        }")
//        println("Root row = ${simplex.M.rows[simplex.M.nRows-1]}")
//        println("Root row +- count ${simplex.M.rows[simplex.M.nRows-1].count {it.value > 0}} ${simplex.M.rows[simplex.M.nRows-1].count {it.value < 0}}")
//        val rootActRows = simplex.M.columns.takeLast(nAgents).map { it.keys.first() }.toSet()
//        println("Count of columns with no root acts ${simplex.M.columns.count {
//            it.all { !rootActRows.contains(it.key) }
//        }}")
//
//        println()
//        println("Starting MCMC chain")
//        for(p in 1..40) {
//            println()
//            val oldSolution = HashIntVector(simplex.X)
//            val pivotTime = measureTimeMillis {
////                simplex.pivotPerturb()
////                simplex.columnPerturb()
////                simplex.hybridPerturb()
//                simplex.mcmcPerturb()
//            }
//            println("Found solution in $pivotTime ms")
//            println("new solution = ${simplex.X}")
//            println("Change is ${simplex.X - oldSolution}")
//            assert(simplex.M*simplex.X == simplex.B)
//            assert(simplex.X.isPositive())
//        }
//
//        println("Done!")
//    }


//    @Test
//    fun findInitialSolutionTest() {
//        val gridSize = 32
//        val timesteps = 6
//        val nAgents = 100
//        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
//        val observations = HashIntVector()
//        val startPos = Array(nAgents) {
//            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
//        }
//        for(agent in startPos.indices) {
//            val xPos = startPos[agent].first
//            val yPos = startPos[agent].second
//            observations[xPos*gridSize + yPos] = -1
//            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
//        }
//
//        val simplex = Simplex(abmMatrix, observations)
//
//
//        println("Starting Hermite transformation")
//        val hermiteTime = measureTimeMillis {
////            simplex.pivotSourcePolyTree()
//            simplex.pivotRootedSourceTree()
//        }
//
//        println("Done in $hermiteTime ms")
//
//        println("pivoted rows = ${simplex.pivotPoints.size}")
//        println("Total M size = ${simplex.M.nRows} x ${simplex.M.nCols}")
//        println("Base solution = ${simplex.X}")
////        println("Observation error = ${simplex.observationNorm()}")
//        assert(simplex.M*simplex.X == simplex.B)
//
//        for(i in 0 until simplex.M.nRows) {
//            assert(simplex.X[simplex.pivotPoints[i]] == simplex.B[i])
//        }
//
//        if(simplex.X.isPositive()) {
//            println("Positive solution already!")
//        } else {
//            println("Pivoting out negatives...")
//            simplex.pivotOutNegatives()
//        }
//
//        assert(simplex.M*simplex.X == simplex.B)
//        assert(simplex.X.isPositive())
//
//    }

//    fun constructInitialSolution(S: OldSimplex, agentPos: Array<Pair<Int,Int>>, gridSize: Int, timesteps: Int): Map<Int,Int> {
//        val pivots = HashMap<Int,Int>()
//        for(agent in agentPos) {
//            for(t in 0 until timesteps step 2) {
//                val startRow = gridSize*gridSize*t + agent.first * gridSize + agent.second
//                val midRow = startRow + gridSize * gridSize + 1
//                val endRow = midRow + gridSize * gridSize - 1
//                val left =
//                    S.M.rows[startRow]
//                        .keys
//                        .find { S.M[midRow, it] == 1 && S.M.columns[it].sparseSize == 2 }
//                        ?: throw(RuntimeException("Action not found"))
//                pivots[midRow] = left
//                val right =
//                    S.M.rows[midRow]
//                        .keys
//                        .find { S.M[endRow, it] == 1 && S.M.columns[it].sparseSize == 2}
//                        ?: throw(RuntimeException("Action not found"))
//                pivots[endRow] = right
//            }
//        }
//        return pivots
//    }
}