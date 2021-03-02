package lib.sparseVector

import lib.abstractAlgebra.IntOperators

inline class IntMapVector(override val nonZeroEntries: MutableMap<Int,Int> = HashMap(4)):
    MutableSparseVector<Int>, IntOperators {
    override fun new() = IntMapVector()
//    override fun toString() = nonZeroEntries.toString()
}