package lib.vector

import lib.abstractAlgebra.Field
import lib.abstractAlgebra.FieldElement
import lib.abstractAlgebra.FieldOperators


interface SparseVector<T: Any>: FieldOperators<T> {
//    val field: Field<T>
    val nonZeroEntries: Map<Int,T>
    fun new(): MutableSparseVector<T>  // place to put results of arithmetic

    // Default implementations

    operator fun get(index: Int): T= nonZeroEntries[index]?:zero
    operator fun plus(other: SparseVector<T>): SparseVector<T> = this.mapTo(new(), other, { it }, { it }, {a,b -> a + b})
    operator fun minus(other: SparseVector<T>): SparseVector<T> = this.mapTo(new(), other, { it }, { -it }, {a,b -> a - b})
    operator fun times(scalar: T): SparseVector<T> = this.mapTo(new()) { it *scalar }
    operator fun unaryMinus(): SparseVector<T> = this.mapTo(new()) { -it }
    fun dotProduct(other: SparseVector<T>): T {
        var sum: T = zero

        for(entry in nonZeroEntries) {
            other.nonZeroEntries[entry.key]?.also { otherValue ->
                sum += entry.value * otherValue
            }
        }
        return sum
    }

    fun mapTo(result: MutableSparseVector<T>, unaryOp: (T) -> T): MutableSparseVector<T> {
        for(entry in nonZeroEntries) {
            result[entry.key] = unaryOp(entry.value)
        }
        return result
    }

    fun mapTo(result: MutableSparseVector<T>, other: SparseVector<T>, lhsOnlyOp: (T)->T, rhsOnlyOp: (T)->T, binaryOp: (T, T) -> T): MutableSparseVector<T> {
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

}


