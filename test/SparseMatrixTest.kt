import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.IntOperators
import lib.sparseMatrix.GridMatrix
import lib.sparseMatrix.times
import lib.sparseVector.MutableMapVector
import lib.sparseVector.asMutableDoubleVector
import lib.sparseVector.asMutableIntVector
import org.junit.Test

class SparseMatrixTest {

    @Test
    fun multiplicationTest() {
        val M = GridMatrix(3,3) { HashMap<Int,Int>().asMutableIntVector() }
        M[0,0] = 2
        M[1,0] = 3
        M[0,2] = 4
        M[2,2] = 1

        val V = MutableMapVector(IntOperators)
        V[0] = 2
        V[2] = -1
        val MV = M * V
        println(MV)
        assert(MV[0] == 0)
        assert(MV[1] == 6)
        assert(MV[2] == -1)
    }

    @Test
    fun multiplicationTestDouble() {
        println(-1.0*0.0 == DoubleOperators.zero)
        val M = GridMatrix(3,3) { HashMap<Int,Double>().asMutableDoubleVector() }
        M[0,0] = 2.0
        M[1,0] = 3.0
        M[0,2] = 4.0
        M[2,2] = 1.0

        val V = MutableMapVector(DoubleOperators)
        V[0] = 2.0
        V[2] = -1.0
        val MV = M * V
        println(MV)
        assert(MV[0] == 0.0)
        assert(MV[1] == 6.0)
        assert(MV[2] == -1.0)
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
