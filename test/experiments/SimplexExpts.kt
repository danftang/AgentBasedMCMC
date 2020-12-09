package experiments

import ABMMatrices.twoDabmMatrix
import Simplex
import SimplexMCMC
import lib.abstractAlgebra.*
import lib.sparseMatrix.GridMapMatrix
import lib.sparseMatrix.IPsolve
import lib.sparseMatrix.mapNonZeroEntriesTo
import lib.vector.*
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.math.absoluteValue
import kotlin.math.pow
import kotlin.random.Random

// Experiments with pivoting the ABM matrix, a-la Simplex algorithm
class SimplexExpts {

    val gridSize = 3
    val timesteps = 2
    val nAgents = 4
    val abmMatrix = twoDabmMatrix(gridSize, timesteps)
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
        simplex.minimise()
        println()
        println(simplex.X())
        assert(simplex.X().nonZeroEntries == mapOf(0 to 10.0, 1 to 30.0, 3 to 20.0))
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
        simplex.minimise()
        println()
        println(simplex.X())
        assert(simplex.X().nonZeroEntries == mapOf(0 to Fraction(10.0), 1 to Fraction(30.0), 3 to Fraction(20.0)))
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


    @Test
    fun twoDABM() {

        val fieldOperators = FractionOperators // DoubleOperators
        // val doubleToField: Double.()->Fraction = { Fraction(this) }
        val intToField: Int.()->Fraction = { Fraction(this,1) }

        val eventProbs = (0 until abmMatrix.nCols).associateWith { fieldOperators.one }.asMapVector(fieldOperators)

        val abmGridMatrix = abmMatrix.mapNonZeroEntriesTo(
            GridMapMatrix(fieldOperators, abmMatrix.nRows, abmMatrix.nCols)) {
            it.intToField()
        }


        val fieldObservations = observations.mapNonZeroEntriesTo(MutableMapVector(fieldOperators)) {
            it.intToField()
        }

//        val initialSolution = abmGridMatrix
//            .IPsolve(fieldObservations, emptyMap<Int,Fraction>().asMapVector(fieldOperators), "==")
//            .toMapVector()
//            .mapNonZeroEntriesTo(MutableMapVector(fieldOperators)) { it.doubleToField() }

//        println(doubleAbmMatrix)
//        println(initialSolution)
//        println(doubleAbmMatrix * initialSolution)
//        println(observations)

//        println("initial solution is ${initialSolution}")

        val simplex = Simplex(abmGridMatrix, fieldObservations, eventProbs)

//        for(i in simplex.basicColsByRow.indices) {
//            val j = simplex.basicColsByRow[i]
//            assert(j >= 0 && j < simplex.bColumn)
//            assert(simplex.M[i,j] != fieldOperators.zero)
//        }

        var nIntegerSolutions = 0
        var nPivots = 0
        for(i in 1..100) {
            println("Sampling $i")
            var ndPivots = simplex.allNonDegeneratePivots()
//            while(ndPivots.isEmpty()) {
            for(j in 1..100) {
                val pivot = simplex.allDegeneratePivots().random()
                simplex.pivot(pivot.row, pivot.col)
                ndPivots = simplex.allNonDegeneratePivots()
            }
            val ndPivot = ndPivots.random()
            simplex.pivot(ndPivot.row, ndPivot.col)
            val x = simplex.X()
            println("Solution is $x")
            if(simplex.isAtIntegerSolution()) {
                println("Solution is Integer")
                nIntegerSolutions++
            } else println("Solution is not integer")

//            if(simplex.isBinary()) println("Binary") else println("Not Binary")
//            if(simplex.isNormalised()) println("Normalised") else println("Not Normalised")
            nPivots++
            println("Fraction of solutions integer = ${nIntegerSolutions*1.0/nPivots}")
//            simplex.setToZeroIfBelow(0.001)
        }


        println("Done!")
    }

//    fun Simplex<Double>.isAtIntegerSolution(): Boolean {
//        val doubleSolution = X()
//        val intSolution = MutableMapVector(IntOperators)
//        doubleSolution.mapNonZeroEntriesTo(intSolution) { it.roundToInt() }
//        return abmMatrix * intSolution == observations
//    }

    fun<T> Simplex<T>.isAtIntegerSolution(): Boolean where T : Comparable<T>, T: Number {
        val intSolution = X().mapNonZeroEntriesTo(MutableMapVector(IntOperators)) { it.roundToInt() }
        return abmMatrix * intSolution == observations
    }


    fun Simplex<Double>.setToZeroIfBelow(smallest: Double) {
        val iter = M.nonZeroEntries.iterator()
        while(iter.hasNext()) {
            val entry = iter.next()
            if(entry.value.absoluteValue < smallest) iter.remove()
        }
    }

}