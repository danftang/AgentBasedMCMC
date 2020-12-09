package lib.vector

import lib.abstractAlgebra.DoubleOperators

inline class DoubleMapVector(override val nonZeroEntries: MutableMap<Int,Double> = HashMap(4)):
    MutableSparseVector<Double>, DoubleOperators {
    override fun new() = DoubleMapVector()
//    override fun toString() = nonZeroEntries.toString()
}