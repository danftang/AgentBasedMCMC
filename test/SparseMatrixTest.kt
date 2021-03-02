import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.IntOperators
import lib.sparseMatrix.GridMapMatrix
import lib.sparseVector.MutableMapVector
import org.junit.Test

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

//    @Test
//    fun IPSolveTest() {
//        val M = GridMapMatrix(DoubleOperators,3,3)
//        M[0,0] = 2.0
//        M[1,0] = 3.0
//        M[0,2] = 4.0
//        M[2,2] = 1.0
//
//        println(M)
//        for(row in M.rows) {
//            for(entry in row.nonZeroEntries) {
//                print("${entry.key} -> ${entry.value}, ")
//            }
//            println()
//        }
//        println()
//        for(col in M.columns) {
//            for(entry in col.nonZeroEntries) {
//                print("${entry.key} -> ${entry.value}, ")
//            }
//            println()
//        }
//
//        println()
//        for(entry in M.nonZeroEntries) {
//            print("${entry.row}:${entry.col} -> ${entry.value}, ")
//        }
//
//
//        val V = MutableMapVector(DoubleOperators)
//        V[0] = 2.0
//        V[2] = 1.0
//        val B = M * V
//        println(B)
//
//        val X = M.IPsolve(B, emptyMap<Int,Double>().asDoubleVector(), "==")
//        println(X.toList())
//    }
}
