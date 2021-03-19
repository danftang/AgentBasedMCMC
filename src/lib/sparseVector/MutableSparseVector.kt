package lib.sparseVector

interface MutableSparseVector<T>: SparseVector<T> {
    // override val nonZeroEntries: MutableMap<Int,T>

    operator fun set(index: Int, value: T)

    fun setToZero()

    fun mapAssign(elementTransform: (T) -> T)

    fun mapAssignWithIndex(transform: (Int, T) -> T)


    // Default implementations
    /////////////////////////////////////
    operator fun timesAssign(multiplier: T) {
        if (multiplier == one) return
        if (multiplier == zero) setToZero()
        mapAssign { it * multiplier }
    }

    operator fun divAssign(denominator: T) {
        if (denominator == one) return
        mapAssign { it / denominator }
    }

    operator fun plusAssign(other: SparseVector<T>) {
        other.nonZeroEntries.forEach { (index, value) ->
            set(index, get(index) + value)
        }
    }

    operator fun minusAssign(other: SparseVector<T>) {
        other.nonZeroEntries.forEach { (index, value) ->
            set(index, get(index) - value)
        }
    }

    fun weightedPlusAssign(other: SparseVector<T>, weight: T) {
        other.nonZeroEntries.forEach { (index, value) ->
            set(index, get(index) + weight * value)
        }
    }
}