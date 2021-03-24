package lib.sparseVector

import lib.abstractAlgebra.DoubleOperators
import lib.abstractAlgebra.FieldOperators
import lib.abstractAlgebra.FractionOperators
import lib.abstractAlgebra.IntOperators
import org.apache.commons.math3.fraction.Fraction

class MutableMapVector<T>(field: FieldOperators<T>, override val nonZeroEntries: MutableMap<Int, T> = HashMap()):
    IMutableMapVector<T>, FieldOperators<T> by field.operators {
    override fun new() = MutableMapVector(operators, HashMap())
    override fun toString() = nonZeroEntries.toString()
    override fun equals(other: Any?) = isEqualTo(other)
}

interface IMutableMapVector<T> : MutableSparseVector<T> {

    override val nonZeroEntries: MutableMap<Int,T>

//    constructor(copyFrom: SparseVector<T>): this(copyFrom.operators, HashMap(copyFrom.nonZeroEntries))

//    override fun toString() = nonZeroEntries.toString()
//    override fun equals(other: Any?) = isEqualTo(other)
    override fun set(index: Int, value: T) { if(value.isZero()) nonZeroEntries.remove(index) else nonZeroEntries[index] = value }
    override fun setToZero() { nonZeroEntries.clear() }


    override fun mapAssign(elementTransform: (T) -> T) {
        val iter = nonZeroEntries.iterator()
        while (iter.hasNext()) {
            val entry = iter.next()
            val newVal = elementTransform(entry.value)
            if (newVal.isZero()) iter.remove() else entry.setValue(newVal)
        }
    }


    override fun mapAssignWithIndex(transform: (Int, T) -> T) {
        val iter = nonZeroEntries.iterator()
        while (iter.hasNext()) {
            val entry = iter.next()
            val newVal = transform(entry.key, entry.value)
            if (newVal.isZero()) iter.remove() else entry.setValue(newVal)
        }
    }


    fun mapAssign(index: Int, remappingFunction: (T)->T) {
        nonZeroEntries.compute(index) { _, oldValue ->
            val newVal = if(oldValue == null) remappingFunction(zero) else remappingFunction(oldValue)
            if(newVal.isZero()) null else newVal
        }
    }


//    override operator fun timesAssign(multiplier: T) {
//        if (multiplier == one) return
//        if (multiplier == zero) nonZeroEntries.clear()
//        mapAssign { it * multiplier }
//    }
//
//    override operator fun divAssign(denominator: T) {
//        if(denominator == one) return
//        mapAssign { it / denominator }
//    }
//
//    override operator fun plusAssign(other: SparseVector<T>) {
//        mapAssign(other, { it }, { it }, {a,b -> a + b})
//    }
//
//    override operator fun minusAssign(other: SparseVector<T>) {
//        mapAssign(other, { it }, { -it }, {a,b -> a - b})
//    }
//
//    override fun weightedPlusAssign(other: SparseVector<T>, weight: T) {
//        mapAssign(other, { it }, { it * weight }, { a, b -> a + b * weight })
//    }


//    fun mapAssign(other: SparseVector<T>, lhsOnlyOp: (T) -> T, rhsOnlyOp: (T) -> T, binaryOp: (T, T) -> T) {
//        val toRemove = ArrayList<Int>()
//        for (entry in nonZeroEntries) {
//            val otherVal = other.nonZeroEntries[entry.key]
//            val newVal = if (otherVal != null) binaryOp(entry.value, otherVal) else lhsOnlyOp(entry.value)
//            if (newVal.isZero()) toRemove.add(entry.key) else entry.setValue(newVal)
//        }
//        for (entry in other.nonZeroEntries) {
//            nonZeroEntries.compute(entry.key) { _, oldValue ->
//                if (oldValue != null) oldValue else {
//                    val newVal = rhsOnlyOp(entry.value)
//                    if (newVal.isZero()) null else newVal
//                }
//            }
//        }
//        toRemove.forEach { nonZeroEntries.remove(it) }
//    }
//

}


inline class MutableDoubleMapVector(override val nonZeroEntries: MutableMap<Int, Double>): IMutableMapVector<Double>, DoubleOperators {
    override fun new() = MutableDoubleMapVector(HashMap())
}
inline fun MutableMap<Int,Double>.asMutableDoubleVector(): MutableDoubleMapVector = MutableDoubleMapVector(this)

inline class MutableFractionMapVector(override val nonZeroEntries: MutableMap<Int, Fraction>): IMutableMapVector<Fraction>,
    FractionOperators {
    override fun new() = MutableFractionMapVector(HashMap())
}
inline fun MutableMap<Int, Fraction>.asMutableFractionVector(): MutableFractionMapVector = MutableFractionMapVector(this)


inline fun<T> MutableMap<Int,T>.asMutableVector(operators: FieldOperators<T>): MutableMapVector<T> =
    MutableMapVector(operators, this)


inline fun MutableMap<Int,Int>.asMutableIntVector(): MutableMapVector<Int> =
    MutableMapVector(IntOperators, this)
