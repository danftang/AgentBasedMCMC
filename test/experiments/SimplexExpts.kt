package experiments

import MutableConstraint
import GridMapSimplex
import lib.abstractAlgebra.*
import lib.sparseVector.*
import org.apache.commons.math3.util.OpenIntToDoubleHashMap
import org.junit.Test
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
    fun anIntroToLinearProgrammingAndGameTheoryPage89() {
        val constraints = listOf(
            MutableConstraint(hashMapOf(0 to -2.0, 2 to 6.0, 3 to 2.0, 5 to -3.0, 6 to 1.0),"==", 20.0),
            MutableConstraint(hashMapOf(0 to -4.0, 1 to 1.0, 2 to 7.0, 3 to 1.0, 5 to -1.0),"==", 10.0),
            MutableConstraint(hashMapOf(2 to -5.0, 3 to 3.0, 4 to 1.0, 5 to -1.0),"==", 60.0)
        )

//        val objective = hashMapOf(2 to 13.0, 3 to -6.0, 5 to 2.0).asMutableDoubleVector()
        val objective = ApacheMapVector(2 to 13.0, 3 to -6.0, 5 to 2.0)
        val initialSolution = hashMapOf(6 to 20.0, 1 to 10.0, 4 to 60.0).asMutableDoubleVector()
        val simplex = GridMapSimplex(constraints, objective, initialSolution)
        println(simplex.M)
        simplex.greedyMinimise()
        println()
        val answer = simplex.X()
        println(answer)
        assert(
            answer.nonZeroEntries == mapOf(0 to 10.0, 1 to 30.0, 3 to 20.0) ||
            answer.nonZeroEntries == mapOf(0 to 1.25, 3 to 22.5, 5 to 7.5)
        )
    }

    @Test
    fun testPivotOutNegatives() {
        val constraints = listOf(
            MutableConstraint(hashMapOf(
                0 to 1.0
            ),"<=", 2.0),
            MutableConstraint(hashMapOf(
                1 to 1.0
            ),"<=", 2.0),
            MutableConstraint(hashMapOf(
                0 to 1.0,
                1 to 1.0
            ), ">=", 1.0)
        )
        val objective = hashMapOf(
            0 to 1.0,
            1 to 1.0
        ).asDoubleVector()
        val simplex = GridMapSimplex(constraints, objective)
//        simplex.pivotToInitialSolutionWithoutORTools()
        println(simplex.M)
        simplex.greedyMinimise()
        println("Solution is")
        println(simplex.X())
        simplex.X(true).nonZeroEntries.forEach { assert(it.value > 0.0) }
    }



//    @Test
//    fun fractionalMatrix() {
//        val coeffs = GridMapMatrix(FractionOperators,4,8)
//        coeffs[0,0] = Fraction(-2.0)
//        coeffs[0,2] = Fraction(6.0)
//        coeffs[0,3] = Fraction(2.0)
//        coeffs[0,5] = Fraction(-3.0)
//        coeffs[0,6] = Fraction(1.0)
//
//        coeffs[1,0] = Fraction(-4.0)
//        coeffs[1,1] = Fraction(1.0)
//        coeffs[1,2] = Fraction(7.0)
//        coeffs[1,3] = Fraction(1.0)
//        coeffs[1,5] = Fraction(-1.0)
//
//        coeffs[2,2] = Fraction(-5.0)
//        coeffs[2,3] = Fraction(3.0)
//        coeffs[2,4] = Fraction(1.0)
//        coeffs[2,5] = Fraction(-1.0)
//
//        // constants
//        coeffs[0,7] = Fraction(20.0)
//        coeffs[1,7] = Fraction(10.0)
//        coeffs[2,7] = Fraction(60.0)
//
//        // objective
//        coeffs[3,2] = Fraction(13.0)
//        coeffs[3,3] = Fraction(-6.0)
//        coeffs[3,5] = Fraction(2.0)
//
//        val initialPivots = intArrayOf(6,1,4)
//        val simplex = Simplex(coeffs, initialPivots)
//        println(simplex.M)
//        simplex.greedyMinimise()
//        println()
//        println(simplex.X())
//        assert(simplex.X().nonZeroEntries == mapOf(0 to Fraction(10.0), 1 to Fraction(30.0), 3 to Fraction(20.0)))
//    }



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




}