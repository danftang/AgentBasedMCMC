package lib.sparseVector

import lib.abstractAlgebra.*
import org.apache.commons.math3.FieldElement

abstract class MapVector<T>(override val nonZeroEntries: Map<Int,T>) : SparseVector<T> {

    override fun equals(other: Any?) = isEqualTo(other)
    override fun toString() = nonZeroEntries.toString()
}

fun DoubleArray.toMapVector(): MapVector<Double> = this.asSequence()
    .withIndex()
    .filter { it.value != 0.0 }
    .associate { Pair(it.index, it.value) }
    .asDoubleVector()


fun<T: FieldElement<T>> Array<T>.toMapVector(field: FieldOperators<T>) = this.asSequence()
    .withIndex()
    .filter { it.value != field.zero }
    .associate { Pair(it.index, it.value) }
    .asVector(field)


inline fun Map<Int,Double>.asDoubleVector(): MapVector<Double> = object: MapVector<Double>(this@asDoubleVector), DoubleOperators {
    override fun new() = DoubleMapVector()
}

inline fun Map<Int,Int>.asIntVector(): MapVector<Int> = object: MapVector<Int>(this@asIntVector), IntOperators {
    override fun new() = IntMapVector()
}

inline fun<T> Map<Int,T>.asVector(operators: FieldOperators<T>): MapVector<T> = object: MapVector<T>(this@asVector),
        FieldOperators<T> by operators {
    override fun new()= MutableMapVector(operators)
}

//inline fun<T: FieldElement<T>> Map<Int,T>.asVector(apacheField: Field<T>): MapVector<T> = asVector(FieldElementOperators(apacheField))

//inline fun Map<Int,DoubleFieldElement>.asDoubleFieldMapVector() = object: MapVector<DoubleFieldElement>(this@asDoubleFieldMapVector) {
//    override fun new(): MutableSparseVector<DoubleFieldElement> {
//        TODO("Not yet implemented")
//    }
//
//    override val zero: DoubleFieldElement
//        get() = DoubleFieldElement(0.0)
//    override val one: DoubleFieldElement
//        get() = DoubleFieldElement(1.0)
//
//}


