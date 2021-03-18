package lib.sparseVector

import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.FieldOperators
import lib.abstractAlgebra.FractionOperators
import lib.abstractAlgebra.IntOperators
import org.apache.commons.math3.fraction.Fraction

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

inline class MutableDoubleMapVector(override val nonZeroEntries: MutableMap<Int, Double>): MutableSparseVector<Double>, DoubleOperators {
    override fun new() = MutableDoubleMapVector(HashMap())
}
inline fun MutableMap<Int,Double>.asMutableDoubleVector(): MutableDoubleMapVector = MutableDoubleMapVector(this)

inline class MutableFractionMapVector(override val nonZeroEntries: MutableMap<Int, Fraction>): MutableSparseVector<Fraction>,
    FractionOperators {
    override fun new() = MutableFractionMapVector(HashMap())
}
inline fun MutableMap<Int, Fraction>.asMutableFractionVector(): MutableFractionMapVector = MutableFractionMapVector(this)


inline fun<T> MutableMap<Int,T>.asMutableVector(operators: FieldOperators<T>): MutableMapVector<T> =
    MutableMapVector(operators, this)


inline fun MutableMap<Int,Int>.asMutableIntVector(): MutableMapVector<Int> =
    MutableMapVector(IntOperators, this)
