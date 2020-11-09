package treeMCMC

import lib.sparseMatrix.HashIntVector
import lib.sparseMatrix.HashColIntMatrix
import lib.sparseMatrix.IntVector
import java.lang.RuntimeException

class PolyhedralCompletionOracle(val M: HashColIntMatrix, val B: HashIntVector): CompletionOracle {

    override fun completion(prefix: List<Int>, breakpointValue: Int): IntVector? {
        val Mprefix = M.subMatrixView(0, prefix.size+1)
        val Mcompletion = M.subMatrixView(prefix.size+1, M.nCols)
        val newValue = IntVector(M.nCols) {  i ->
            if(i < prefix.size) prefix[i] else 0
        }
        newValue[prefix.size] = breakpointValue

        val target = B - Mprefix*IntVector(newValue.subList(0, prefix.size+1)).toSparseIntVector()
        return try {
            val completion = Mcompletion.IPsolve(target, DoubleArray(Mcompletion.nCols) { 0.0 }.asList(), "=")
            for (i in prefix.size + 1 until newValue.size) {
                newValue[i] = completion[i]
            }
            newValue
        } catch(exception: RuntimeException) {
            null
        }
    }


    override fun upperBounds(solution: IntArray, prefixLength: Int): IntArray {
        TODO("Not yet implemented")
    }

}