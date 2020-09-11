package treeMCMC

import lib.SparseColIntMatrix
import java.lang.RuntimeException

class PolyhedralCompletionOracle(val M: SparseColIntMatrix, val B: SparseColIntMatrix.SparseIntColumn): CompletionOracle {

    override fun completion(prefix: List<Int>, breakpointValue: Int): IntArray? {
        val Mprefix = M.subMartixView(0, prefix.size+1)
        val Mcompletion = M.subMartixView(prefix.size+1, M.nCols)
        val newValue = IntArray(M.nCols) {  i ->
            if(i < prefix.size) prefix[i] else 0
        }
        newValue[prefix.size] = breakpointValue

        val target = B - Mprefix*newValue.asList().subList(0, prefix.size+1)
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