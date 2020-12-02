package lib.vector

import lib.abstractAlgebra.FieldElement
import org.apache.commons.math3.linear.OpenMapRealVector
import org.apache.commons.math3.linear.SparseFieldVector

interface MutableSparseVector<T: Any>: SparseVector<T> {
    override val nonZeroEntries: MutableMap<Int,T>

    operator fun set(index: Int, value: T){
        if(value != zero) nonZeroEntries[index] = value else nonZeroEntries.remove(index)
    }

    fun mapAssign(index: Int, remappingFunction: (T)->T) {
        nonZeroEntries.compute(index) { _, oldValue ->
            val newVal = if(oldValue == null) remappingFunction(zero) else remappingFunction(oldValue)
            if(newVal == zero) null else newVal
        }
    }

    operator fun timesAssign(multiplier: T) {
        if (multiplier == one) return
        if (multiplier == zero) nonZeroEntries.clear()
        mapAssign { it * multiplier }
    }

    operator fun divAssign(denominator: T) {
        if(denominator == one) return
        mapAssign { it / denominator }
    }

    operator fun plusAssign(other: SparseVector<T>) {
        mapAssign(other, { it }, { it }, {a,b -> a + b})
    }

    operator fun minusAssign(other: SparseVector<T>) {
        mapAssign(other, { it }, { -it }, {a,b -> a - b})
    }

    fun weightedPlusAssign(other: SparseVector<T>, weight: T) {
        mapAssign(other, { it }, { it * weight }, { a,b -> a + b * weight } )
    }

}

inline fun<T: Any> MutableSparseVector<T>.mapAssign(unaryOp: (T)->T) {
    val iter = nonZeroEntries.iterator()
    while(iter.hasNext()) {
        val entry = iter.next()
        val newVal = unaryOp(entry.value)
        if(newVal != zero) entry.setValue(newVal) else iter.remove()
    }
}

inline fun<T: Any> MutableSparseVector<T>.mapAssign(other: SparseVector<T>, lhsOnlyOp: (T)->T, crossinline rhsOnlyOp: (T)->T, binaryOp: (T, T)->T) {
    val iter = nonZeroEntries.iterator()
    while(iter.hasNext()) {
        val entry = iter.next()
        val otherVal = other.nonZeroEntries[entry.key]
        val newVal = if(otherVal != null) binaryOp(entry.value, otherVal) else lhsOnlyOp(entry.value)
        if(newVal != zero) entry.setValue(newVal) else iter.remove()
    }
    for(entry in other.nonZeroEntries) {
        nonZeroEntries.compute(entry.key) { key, oldValue ->
            if(oldValue != null) oldValue else {
                val newVal = rhsOnlyOp(entry.value)
                if(newVal == zero) null else newVal
            }
        }
    }
}
