import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.FractionOperators
import lib.abstractAlgebra.IntOperators
import lib.collections.GridMap
import lib.plus
import lib.sparseMatrix.GridMapMatrix
import lib.sparseVector.MutableMapVector
import org.apache.commons.math3.fraction.Fraction
import org.junit.Test
import kotlin.random.Random
import kotlin.system.measureTimeMillis

class SparseMatrixTest {

//    init {
//        System.loadLibrary("jniortools")
//    }

    @Test
    fun multiplicationTest() {
        val M = GridMapMatrix(IntOperators,3,3)
        M[0,0] = 2
        M[1,0] = 3
        M[0,2] = 4
        M[2,2] = 1

        val V = MutableMapVector(IntOperators)
        V[0] = 2
        V[2] = -1
        println(M * V)
    }

    @Test
    fun multiplicationTestDouble() {
        println(-1.0*0.0 == DoubleOperators.zero)
        val M = GridMapMatrix(DoubleOperators,3,3)
        M[0,0] = 2.0
        M[1,0] = 3.0
        M[0,2] = 4.0
        M[2,2] = 1.0

        val V = MutableMapVector(DoubleOperators)
        V[0] = 2.0
        V[2] = -1.0
        println(M * V)
    }


    @Test
    fun gridMapMatrixStressTest() {
        println("Setting up matrix")
        val M = GridMapMatrix(FractionOperators, 1000,5000)
        for(i in 0 until M.nRows) {
            for(j in 0 until M.nCols) {
                if(Random.nextDouble() < 0.5) {
                    M[i, j] = Fraction.ONE
                }
            }
        }
        println("Doing transformations")
        var tot = 0.0
        val tim = measureTimeMillis {
            for (n in 1..50000000) {
                val r = Random.nextInt(M.nRows)
                val c = Random.nextInt(M.nCols)
//            M[r, c] += Fraction(1,2)
                tot += M[r, c].toDouble()
//                tot += (M.gridMap[r,c]?:M.zero).toDouble()
//            M.mapAssign(Random.nextInt(M.nRows), Random.nextInt(M.nCols)) { it + Fraction(1,2) }
            }
        }
        println(tot)
        println(tim/1000.0)
    }

    @Test
    fun gridMapStressTest() {
        println("Setting up matrix")
        val M = GridMap<Fraction>( 1000,5000)
        for(i in 0 until M.nRows) {
            for(j in 0 until M.nCols) {
                if(Random.nextDouble() < 0.5) {
                    M[i, j] = Fraction.ONE
                }
            }
        }
        println("Doing transformations")
        var tot = 0.0
        val tim = measureTimeMillis {
            for (n in 1..50000000) {
                val r = Random.nextInt(M.nRows)
                val c = Random.nextInt(M.nCols)
                tot += (M[r, c]?: Fraction.ZERO).toDouble()
            }
        }
        println(tot)
        println(tim/1000.0)
    }

}
