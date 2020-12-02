package treeMCMC

import lib.sparseIntMatrix.HashIntVector
import lib.sparseIntMatrix.HashColIntMatrix
import lib.sparseIntMatrix.IntArrayVector
import java.lang.RuntimeException

class PolyhedralCompletionOracle(val M: HashColIntMatrix, val B: HashIntVector): CompletionOracle {

    override fun completion(prefix: List<Int>, breakpointValue: Int): IntArrayVector? {
        val Mprefix = M.subMatrixView(0, prefix.size+1)
        val Mcompletion = M.subMatrixView(prefix.size+1, M.nCols)
        val newValue = IntArrayVector(M.nCols) { i ->
            if(i < prefix.size) prefix[i] else 0
        }
        newValue[prefix.size] = breakpointValue

        val target = B - Mprefix* IntArrayVector(newValue.subList(0, prefix.size+1)).toSparseIntVector()
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