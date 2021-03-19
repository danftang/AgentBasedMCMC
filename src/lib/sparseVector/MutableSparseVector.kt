package lib.sparseVector

interface MutableSparseVector<T>: SparseVector<T> {
    // override val nonZeroEntries: MutableMap<Int,T>

    fun remove(index: Int)
    operator fun set(index: Int, value: T)
    fun clear()

    fun mapAssign(elementTransform: (T) -> T)

    operator fun timesAssign(multiplier: T) {
        if (multiplier == one) return
        if (multiplier == zero) clear()
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