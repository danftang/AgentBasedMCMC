package experiments

import Constraint
import Simplex
import lib.abstractAlgebra.*
import lib.sparseMatrix.GridMapMatrix
import lib.vector.*
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.math.absoluteValue
import kotlin.random.Random

// Experiments with pivoting the ABM matrix, a-la Simplex algorithm
class SimplexExpts {

    val gridSize = 32
    val timesteps = 4
    val nAgents = 30
    val observations = MutableMapVector(IntOperators)


    init {
        System.loadLibrary("jniortools")

        val startPos = Array(nAgents) {
            Pair(Random.nextInt(gridSize-1), Random.nextInt(gridSize-1))
        }
        for(agent in startPos.indices) {
            val xPos = startPos[agent].first
            val yPos = startPos[agent].second
            observations[xPos*gridSize + yPos] = -1
            observations[gridSize*gridSize*timesteps + xPos*gridSize + yPos] = 1
        }

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
        simplex.greedyMinimise()
        println()
        println(simplex.X())
        assert(simplex.X().nonZeroEntries == mapOf(0 to 10.0, 1 to 30.0, 3 to 20.0))
    }

    @Test
    fun testPivotOutNegatives() {
        val constraints = listOf(
            Constraint(hashMapOf(
                0 to 1.0
            ),"<=", 2.0),
            Constraint(hashMapOf(
                1 to 1.0
            ),"<=", 2.0),
            Constraint(hashMapOf(
                0 to 1.0,
                1 to 1.0
            ), ">=", 1.0)
        )
        val objective = hashMapOf(
            0 to 1.0,
            1 to 1.0
        ).asDoubleVector()
        val simplex = Simplex(constraints, objective)
        simplex.pivotToInitialSolutionWithoutORTools()
        println(simplex.M)
        simplex.greedyMinimise()
        println("Solution is")
        println(simplex.X())
        simplex.X(true).nonZeroEntries.forEach { assert(it.value > 0.0) }
    }



    @Test
    fun fractionalMatrix() {
        val coeffs = GridMapMatrix(FieldElementOperators(Fraction.ZERO.field),4,8)
        coeffs[0,0] = Fraction(-2.0)
        coeffs[0,2] = Fraction(6.0)
        coeffs[0,3] = Fraction(2.0)
        coeffs[0,5] = Fraction(-3.0)
        coeffs[0,6] = Fraction(1.0)

        coeffs[1,0] = Fraction(-4.0)
        coeffs[1,1] = Fraction(1.0)
        coeffs[1,2] = Fraction(7.0)
        coeffs[1,3] = Fraction(1.0)
        coeffs[1,5] = Fraction(-1.0)

        coeffs[2,2] = Fraction(-5.0)
        coeffs[2,3] = Fraction(3.0)
        coeffs[2,4] = Fraction(1.0)
        coeffs[2,5] = Fraction(-1.0)

        // constants
        coeffs[0,7] = Fraction(20.0)
        coeffs[1,7] = Fraction(10.0)
        coeffs[2,7] = Fraction(60.0)

        // objective
        coeffs[3,2] = Fraction(13.0)
        coeffs[3,3] = Fraction(-6.0)
        coeffs[3,5] = Fraction(2.0)

        val initialPivots = intArrayOf(6,1,4)
        val simplex = Simplex(coeffs, initialPivots)
        println(simplex.M)
        simplex.greedyMinimise()
        println()
        println(simplex.X())
        assert(simplex.X().nonZeroEntries == mapOf(0 to Fraction(10.0), 1 to Fraction(30.0), 3 to Fraction(20.0)))
    }



//    fun Simplex<Double>.isAtIntegerSolution(): Boolean {
//        val doubleSolution = X()
//        val intSolution = MutableMapVector(IntOperators)
//        doubleSolution.mapNonZeroEntriesTo(intSolution) { it.roundToInt() }
//        return abmMatrix * intSolution == observations
//    }

//    fun<T> Simplex<T>.isAtIntegerSolution(): Boolean where T : Comparable<T>, T: Number {
//        val intSolution = X().mapNonZeroEntriesTo(MutableMapVector(IntOperators)) { it.roundToInt() }
//        return abmMatrix * intSolution == observations
//    }


    fun Simplex<Double>.setToZeroIfBelow(smallest: Double) {
        M.nonZeroEntries.forEach { entry ->
            if(entry.value.absoluteValue < smallest) entry.setValue(zero)
        }
//        val iter = M.nonZeroEntries.iterator()
//        while(iter.hasNext()) {
//            val entry = iter.next()
//            if(entry.value.absoluteValue < smallest) iter.remove()
//        }
    }


}