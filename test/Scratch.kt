import lib.Gnuplot
import lib.SettableLazy
import lib.abstractAlgebra.FractionOperators
import lib.gnuplot
import lib.plus
import lib.sparseMatrix.GridMapMatrix
import org.apache.commons.math3.fraction.Fraction
import org.apache.commons.math3.linear.OpenMapRealMatrix
import org.apache.commons.math3.linear.OpenMapRealVector
import org.junit.Test
import kotlin.math.min
import kotlin.random.Random
import kotlin.system.measureTimeMillis

class Scratch {


    @Test
    fun stuff() {
        println("Setting up")

        val M = Array(5000) {
            val col = OpenMapRealVector(1000)
            for(i in 0 until 1000) {
                if(Random.nextDouble() < 0.5) {
                    col.setEntry(i, 1.0)
                }
            }
            col
        }

        println("Doing transformations")
        var tot = 0.0
        val tim = measureTimeMillis {
            for (n in 1..50000000) {
                val r = Random.nextInt(1000)
                val c = Random.nextInt(5000)
                M[c].setEntry(r, M[c].getEntry(r) + 0.5)
                tot += M[c].getEntry(r)
            }
        }
        println(tot)
        println(tim)
        matrix(1.0)
    }



    fun matrix(num: Double) {
        println("Setting up")
        val M = ArrayList<HashMap<Int,Double>>()
        for(j in 0 until 5000) {
            M.add(HashMap())
            for(i in 0 until 1000) {
                if(Random.nextDouble() < 0.5) {
                    M[j][i] = num
                }
            }
        }

        println("Doing transformations")
        var tot = 0.0
        val tim = measureTimeMillis {
            for (n in 1..50000000) {
                val r = Random.nextInt(1000)
                val c = Random.nextInt(5000)
                M[c][r] = M[c][r]?:0.0 + num
                tot += M[c][r] ?: 0.0
            }
        }
        println(tot)
        println(tim)
    }


}