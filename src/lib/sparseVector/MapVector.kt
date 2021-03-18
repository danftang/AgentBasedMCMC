package lib.sparseVector

import lib.abstractAlgebra.*
import org.apache.commons.math3.FieldElement
import org.apache.commons.math3.fraction.Fraction

class MapVector<T>(fieldOperators: FieldOperators<T>, override val nonZeroEntries: Map<Int,T>) :
    SparseVector<T>,
    FieldOperators<T> by fieldOperators {

    override fun new()= MutableMapVector(operators)
    override fun equals(other: Any?) = isEqualTo(other)
    override fun toString() = nonZeroEntries.toString()
}

inline class DoubleMapVector(override val nonZeroEntries: Map<Int, Double>): SparseVector<Double>, DoubleOperators {
    override fun new() = MutableDoubleMapVector(HashMap())
}
inline fun Map<Int,Double>.asDoubleVector(): DoubleMapVector = DoubleMapVector(this)

inline class FractionMapVector(override val nonZeroEntries: Map<Int, Fraction>): SparseVector<Fraction>, FractionOperators {
    override fun new() = MutableFractionMapVector(HashMap())
}
inline fun Map<Int, Fraction>.asFractionVector(): FractionMapVector = FractionMapVector(this)


fun DoubleArray.toDoubleMapVector(): MapVector<Double> = this.asSequence()
    .withIndex()
    .filter { it.value != 0.0 }
    .associate { Pair(it.index, it.value) }
    .asVector(DoubleOperators)

fun IntArray.toIntMapVector(): MapVector<Int> = this.asSequence()
    .withIndex()
    .filter { it.value != 0 }
    .associate { Pair(it.index, it.value) }
    .asVector(IntOperators)


fun<T: FieldElement<T>> Array<T>.toMapVector(field: FieldOperators<T>): MapVector<T> = this.asSequence()
    .withIndex()
    .filter { it.value != field.zero }
    .associate { Pair(it.index, it.value) }
    .asVector(field)



inline fun Map<Int,Int>.asIntVector(): MapVector<Int> = MapVector(IntOperators, this)
inline fun<T> Map<Int,T>.asVector(operators: FieldOperators<T>): MapVector<T> = MapVector(operators, this)
inline fun<T> emptyMapVector(operators: FieldOperators<T>): MapVector<T> = MapVector(operators, emptyMap())
inline fun<T> emptySparseVector(operators: FieldOperators<T>): SparseVector<T> = MapVector(operators, emptyMap())

