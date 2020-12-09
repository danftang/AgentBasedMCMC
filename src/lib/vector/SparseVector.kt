package lib.vector

import lib.abstractAlgebra.AlgebraicOperators
import lib.abstractAlgebra.FieldOperators
import java.lang.RuntimeException
import java.lang.StringBuilder


interface SparseVector<T>: FieldOperators<T> {
//    val field: Field<T>
    val nonZeroEntries: Map<Int,T>
    fun new(): MutableSparseVector<T>  // place to put results of arithmetic

    // Default implementations

    operator fun get(index: Int): T= nonZeroEntries[index]?:zero
    operator fun plus(other: SparseVector<T>): SparseVector<T> = this.mapNonZeroEntriesTo(new(), other, { it }, { it }, { a, b -> a + b})
    operator fun minus(other: SparseVector<T>): SparseVector<T> = this.mapNonZeroEntriesTo(new(), other, { it }, { -it }, { a, b -> a - b})
    operator fun times(scalar: T): SparseVector<T> = this.mapNonZeroEntriesTo(new()) { it *scalar }
    operator fun unaryMinus(): SparseVector<T> = this.mapNonZeroEntriesTo(new()) { -it }
    fun dotProduct(other: SparseVector<T>): T {
        var sum: T = zero

        for(entry in nonZeroEntries) {
            other.nonZeroEntries[entry.key]?.also { otherValue ->
                sum += entry.value * otherValue
            }
        }
        return sum
    }

    fun<R> mapNonZeroEntriesTo(result: MutableSparseVector<R>, unaryOp: (T) -> R): MutableSparseVector<R> {
        for(entry in nonZeroEntries) {
            result[entry.key] = unaryOp(entry.value)
        }
        return result
    }

    fun<R> mapNonZeroEntriesTo(result: MutableSparseVector<R>, other: SparseVector<T>, lhsOnlyOp: (T)->R, rhsOnlyOp: (T)->R, binaryOp: (T, T) -> R): MutableSparseVector<R> {
        for(entry in nonZeroEntries) {
            val otherVal = other.nonZeroEntries[entry.key]
            result[entry.key] =  if(otherVal != null) {
                binaryOp(entry.value, otherVal)
            } else {
                lhsOnlyOp(entry.value)
            }
        }
        for(entry in other.nonZeroEntries) {
            if(!nonZeroEntries.containsKey(entry.key)) result[entry.key] = rhsOnlyOp(entry.value)
        }
        return result
    }

    fun isEqualTo(other: Any?): Boolean {
        return if(other is SparseVector<*>) {
            if(operators != other.operators) throw(RuntimeException("Comparing vectors of different type. You probably didn't want to do that!"))
            operators == other.operators && nonZeroEntries == other.nonZeroEntries
        } else if(other is Vector<*>) {
            (0 until other.size).all { other[it] == this[it] }
        } else false
    }
}

inline fun<T> SparseVector<T>.copyTo(destination: MutableSparseVector<T>): MutableSparseVector<T> {
    for(entry in nonZeroEntries) destination[entry.key] = entry.value
    return destination
}

inline fun<T> SparseVector<T>.toMutableSparseVector(): MutableSparseVector<T> {
    return this.copyTo(new())
}



