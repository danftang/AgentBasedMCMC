package experiments

import ABMMatrices.twoDabmMatrix
import Simplex
import lib.SparseColIntMatrix
import org.junit.Test
import java.lang.RuntimeException
import kotlin.math.absoluteValue
import kotlin.random.Random
import kotlin.system.measureTimeMillis

// Experiments with pivoting the ABM matrix, a-la Simplex algorithm
class SimplexExpts {

    init {
        System.loadLibrary("jniortools")
    }

    @Test
    fun oneDABM() {
        val timesteps = 5
        val N = 1
        val abmMatrix = ABMMatrices.oneDNoInteractionMatrix(timesteps, N)
        val observations = SparseColIntMatrix.SparseIntColumn()
        observations[0] = -1
        observations[timesteps*N] = 1
        abmMatrix.addColRight(SparseColIntMatrix.SparseIntColumn(1 to -1, 3 to -1, 4 to 1))

//        println(abmMatrix)
        val simplex = Simplex(abmMatrix, observations)
        println(simplex)
        simplex.pivotPolyTree()
//        simplex.pivotAt(4,3)
//        println(simplex)
//        simplex.pivotAt(3,2)
//        println(simplex)
//        simplex.pivotAt(2,1)
//        println(simplex)
//        simplex.pivotAt(1,0)
        println(simplex)
        println(simplex.findAllPositivePivotCols())
    }

    @Test
    fun twoDABM() {
        val gridSize = 32
        val timesteps = 6
        val nAgents = 100
        val abmMatrix = twoDabmMatrix(gridSize, timesteps)
        val observations = SparseColIntMatrix.SparseIntColumn()
        val startPos = Array(nAgents) {
            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
        }
        for(agent in startPos.indices) {
            val xPos = startPos[agent].first
            val yPos = startPos[agent].second
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }

        val simplex = Simplex(abmMatrix, observations)
        val Xinit = constructInitialSolution(simplex, startPos, gridSize, timesteps)

        println("Starting Hermite transformation")
        val hermiteTime = measureTimeMillis {
            simplex.pivotPolyTree(Xinit)
        }
        println("Done in $hermiteTime ms")

        println("pivoted rows = ${simplex.pivotPoints.size}")
        println("Total rows = ${simplex.nRows}")
        println("Base solution = ${simplex.X}")
        println("Observation error = ${simplex.observationNorm()}")

        assert(simplex.M*simplex.X == simplex.B)
        assert(simplex.X.isPositive())


        println()
        println("Starting MCMC chain")
        for(p in 1..40) {
            println()
            val pivotTime = measureTimeMillis {
                simplex.pivotPerturb()
            }
            println("Found solution in $pivotTime ms")
            println("new solution = ${simplex.X}")
            println("Observation error = ${simplex.observationNorm()}")
            assert(simplex.M*simplex.X == simplex.B)
            assert(simplex.X.isPositive())
        }

        println("Done!")
    }


    fun constructInitialSolution(S: Simplex, agentPos: Array<Pair<Int,Int>>, gridSize: Int, timesteps: Int): Map<Int,Int> {
        val pivots = HashMap<Int,Int>()
        for(agent in agentPos) {
            for(t in 0 until timesteps step 2) {
                val startRow = gridSize*gridSize*t + agent.first * gridSize + agent.second
                val midRow = startRow + gridSize * gridSize + 1
                val endRow = midRow + gridSize * gridSize - 1
                val left =
                    S.MT[startRow]
                        .keys
                        .find { S.M[midRow, it] == 1 && S.M[it].sparseSize == 2 }
                        ?: throw(RuntimeException("Action not found"))
                pivots[midRow] = left
                val right =
                    S.MT[midRow]
                        .keys
                        .find { S.M[endRow, it] == 1 && S.M[it].sparseSize == 2}
                        ?: throw(RuntimeException("Action not found"))
                pivots[endRow] = right
            }
        }
        return pivots
    }
}