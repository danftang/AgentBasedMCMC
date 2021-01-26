package lib.sparseMatrix

import lib.collections.RemovableEntry
import lib.collections.mutableMapByRemovable
import lib.vector.MutableSparseVector
import lib.vector.SparseVector

interface SparseRowMatrix<T: Any>: SparseMatrix<T> {
    val rows: List<MutableSparseVector<T>>

    override val nRows: Int
        get() = rows.size

    override val nonZeroEntries: Iterable<SparseMatrix.Entry<T>>
        get() = rows.withIndex().flatMap { (row, vector) ->
            vector.nonZeroEntries.entries.mutableMapByRemovable { entry ->
                REntry(row, entry) { this.isZero() }
            }
        }

    override fun get(row: Int, col: Int): T {
        return rows[row][col]
    }

    override fun set(row: Int, col: Int, value: T) {
        rows[row][col] = value
    }

    override fun mapAssign(row: Int, col: Int, remappingFunction: (T) -> T) {
        rows[row].mapAssign(col, remappingFunction)
    }


    override operator fun times(X: SparseVector<T>): SparseVector<T> {
        val result = rows.first().new()
        for(i in rows.indices) {
            result[i] = rows[i].dotProduct(X)
        }
        return result
    }

    class REntry<T>(override val row: Int, val entry: RemovableEntry<MutableMap.MutableEntry<Int, T>>, val isZero: T.()->Boolean): SparseMatrix.Entry<T> {
        override val col: Int
            get() = entry.value.key
        override val value: T
            get() = entry.value.value
        override fun setValue(newValue: T): T {
            return if(newValue.isZero()) {
                val oldValue = value
                entry.remove()
                oldValue
            } else {
                entry.value.setValue(newValue)
            }
        }
    }
}