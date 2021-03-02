package lib.sparseVector

import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.FieldOperators
import lib.abstractAlgebra.IntOperators

class MutableMapVector<T>(
    val fieldOperators: FieldOperators<T>,
    override val nonZeroEntries: MutableMap<Int,T> = HashMap())
    : MutableSparseVector<T>, FieldOperators<T> by fieldOperators
{

    constructor(copyFrom: SparseVector<T>): this(copyFrom.operators, HashMap(copyFrom.nonZeroEntries))

    override fun new() = MutableMapVector(fieldOperators)
    override fun toString() = nonZeroEntries.toString()
    override fun equals(other: Any?) = isEqualTo(other)
}

inline fun<T> MutableMap<Int,T>.asMutableVector(operators: FieldOperators<T>): MutableMapVector<T> =
    MutableMapVector(operators, this)

inline fun MutableMap<Int,Double>.asDoubleMutableVector(): MutableMapVector<Double> =
    MutableMapVector(DoubleOperators, this)

inline fun MutableMap<Int,Int>.asIntMutableVector(): MutableMapVector<Int> =
    MutableMapVector(IntOperators, this)
