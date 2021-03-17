package lib.sparseVector

interface MutableSparseVector<T>: SparseVector<T> {
    override val nonZeroEntries: MutableMap<Int,T>

    operator fun set(index: Int, value: T){
        if(value.isZero()) nonZeroEntries.remove(index) else nonZeroEntries[index] = value
    }

    fun mapAssign(index: Int, remappingFunction: (T)->T) {
        nonZeroEntries.compute(index) { _, oldValue ->
            val newVal = if(oldValue == null) remappingFunction(zero) else remappingFunction(oldValue)
            if(newVal.isZero()) null else newVal
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

inline fun<T> MutableSparseVector<T>.mapAssign(unaryOp: (T)->T) {
    val iter = nonZeroEntries.iterator()
    while(iter.hasNext()) {
        val entry = iter.next()
        val newVal = unaryOp(entry.value)
        if(newVal.isZero()) iter.remove() else entry.setValue(newVal)
    }
}

inline fun<T> MutableSparseVector<T>.mapAssign(other: SparseVector<T>, lhsOnlyOp: (T)->T, crossinline rhsOnlyOp: (T)->T, binaryOp: (T, T)->T) {
    val toRemove = ArrayList<Int>()
    for(entry in nonZeroEntries) {
        val otherVal = other.nonZeroEntries[entry.key]
        val newVal = if(otherVal != null) binaryOp(entry.value, otherVal) else lhsOnlyOp(entry.value)
        if(newVal.isZero()) toRemove.add(entry.key) else entry.setValue(newVal)
    }
    for(entry in other.nonZeroEntries) {
        nonZeroEntries.compute(entry.key) { _, oldValue ->
            if(oldValue != null) oldValue else {
                val newVal = rhsOnlyOp(entry.value)
                if(newVal.isZero()) null else newVal
            }
        }
    }
    toRemove.forEach { nonZeroEntries.remove(it) }
}
