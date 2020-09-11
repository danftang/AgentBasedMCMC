package treeMCMC

interface CompletionOracle {
    fun completion(prefix: List<Int>, breakpointValue: Int): IntArray?
    fun upperBounds(solution: IntArray, prefixLength: Int): IntArray
}
