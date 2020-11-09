package treeMCMC

import lib.sparseMatrix.IntVector

interface CompletionOracle {
    fun completion(prefix: List<Int>, breakpointValue: Int): IntVector?
    fun upperBounds(solution: IntArray, prefixLength: Int): IntArray
}
