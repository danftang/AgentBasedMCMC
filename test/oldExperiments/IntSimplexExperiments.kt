package experiments

import ABMMatrices.twoDabmMatrix
import IntSimplex
import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.HashRowColIntMatrix
import lib.sparseIntMatrix.IntArrayVector
import org.junit.Test
import kotlin.math.absoluteValue
import kotlin.random.Random

class IntSimplexExperiments {

    init {
        System.loadLibrary("jniortools")
    }


    @Test
    fun introToLinearProgrammingAndGameTheoryPage89() {
        val coeffs = HashRowColIntMatrix(3,7)
        coeffs[0,0] = -2
        coeffs[0,2] = 6
        coeffs[0,3] = 2
        coeffs[0,5] = -3
        coeffs[0,6] = 1

        coeffs[1,0] = -4
        coeffs[1,1] = 1
        coeffs[1,2] = 7
        coeffs[1,3] = 1
        coeffs[1,5] = -1

        coeffs[2,2] = -5
        coeffs[2,3] = 3
        coeffs[2,4] = 1
        coeffs[2,5] = -1

        val consts = HashIntVector(0 to 20, 1 to 10, 2 to 60)
        val objective = HashIntVector(2 to 13, 3 to -6, 5 to 2)
//        val initialSolution = HashIntVector(1 to 10, 4 to 60, 6 to 20)
        val simplex = IntSimplex(coeffs, consts, objective)
        simplex.minimise()
        println(simplex.X)
    }


    // Show that we can pivot between any two solutions
    // while maintaining positive, integer intermediate
    // solutions
    @Test
    fun pivotBetweenSolutions() {
        // setup the problem
        val gridSize = 8
        val timesteps = 4
        val nAgents = 20
        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
        val observations = HashIntVector()
        val startPos = Array(nAgents) {
            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
        }
        for(agent in startPos.indices) {
            val xPos = startPos[agent].first
            val yPos = startPos[agent].second
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }


        val Xb = abmMatrix.IPsolve(observations, DoubleArray(abmMatrix.nCols) { Random.nextDouble() }.asList(), "==")
        println(HashIntVector(Xb).sortedBy { it.key })

        // OR-tools solve
//        val doubleObjective = DoubleArray(abmMatrix.nCols) { if(Xb[it] == 0) 1.0 else 0.0 }.asList()
//        val minCorner = abmMatrix.IPsolve(observations, doubleObjective, "==")
//        println(HashIntVector(minCorner).sortedBy { it.key })

        val objective = HashIntVector()
        for(j in 0 until abmMatrix.nCols) {
            if(Xb[j] == 0) objective[j] = 1
        }
        val simplex = IntSimplex(abmMatrix, observations, objective)
        simplex.minimise()
        println(simplex.X)

        println("Done!!!")
    }


    @Test
    fun minimise() {
        val simplex = twoDabmSimplex(HashIntVector(IntArrayVector(IntArray(100) { if(Random.nextBoolean()) 1 else 0 })))
        simplex.minimise()
    }

    // Pivot along a random edge and see the distribution of fractional solutions
    @Test
    fun randomPivot() {
        val simplex = twoDabmSimplex()

        for(i in simplex.basicColsByRow.indices) {
            val j = simplex.basicColsByRow[i]
            assert(j >= 0 && j < simplex.bColumn)
            assert(simplex.M[i,j] != 0)
        }

        var nIntegerSolutions = 0
        var nPivots = 0
        for(i in 1..100) {
            var nNullPivots = -1
            do {
                nNullPivots++
                val pivot = simplex.allPositivePivotPoints().random()
                simplex.pivot(pivot.row, pivot.col)
            } while(simplex.B[pivot.row] == 0)
            println("$nNullPivots null pivots")
            if(simplex.isAtIntegerSolution()) {
                println("Solution is Integer")
                nIntegerSolutions++
            } else println("Solution is not integer")
//            if(simplex.isBinary()) println("Binary") else println("Not Binary")
//            if(simplex.isNormalised()) println("Normalised") else println("Not Normalised")
            nPivots++
            println("Fraction of solutions integer = ${nIntegerSolutions*1.0/nPivots}")
        }
    }

    fun twoDabmSimplex(objective: HashIntVector = HashIntVector()): IntSimplex {
        // setup the problem
        val gridSize = 8
        val timesteps = 6
        val nAgents = 20
        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
        val observations = HashIntVector()
        val startPos = Array(nAgents) {
            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
        }
        for(agent in startPos.indices) {
            val xPos = startPos[agent].first
            val yPos = startPos[agent].second
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }
        return IntSimplex(abmMatrix, observations, objective)
    }

    fun IntSimplex.isBinary(): Boolean {
        return M.entries.all { entry ->
            entry.value.absoluteValue == 1 || entry.col == bColumn || entry.row == objectiveRow
        }
    }

    fun IntSimplex.isNormalised(): Boolean {
        return basicColsByRow.withIndex().all { entry ->
            M[entry.index, entry.value] == 1
        }
    }

}