package treeMCMC

import lib.sparseIntMatrix.IntArrayVector

interface CompletionOracle {
    fun completion(prefix: List<Int>, breakpointValue: Int): IntArrayVector?
    fun upperBounds(solution: IntArray, prefixLength: Int): IntArray
}
