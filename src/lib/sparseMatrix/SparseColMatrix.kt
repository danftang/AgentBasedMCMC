package lib.sparseMatrix

import lib.collections.RemovableEntry
import lib.collections.mutableMapByRemovable
import lib.vector.MutableSparseVector
import lib.vector.SparseVector

interface SparseColMatrix<T: Any>: SparseMatrix<T> {
    val columns: List<MutableSparseVector<T>>

    override val nCols: Int
        get() = columns.size

    override val nonZeroEntries: Iterable<SparseMatrix.Entry<T>>
        get() = columns.withIndex().flatMap { (j, col) ->
            col.nonZeroEntries.entries.mutableMapByRemovable { entry ->
                CEntry(j, entry) { this.isZero() }
            }
        }

    override fun get(row: Int, col: Int): T {
        return columns[col][row]
    }

    override fun set(row: Int, col: Int, value: T) {
        columns[col][row] = value
    }

    override fun mapAssign(row: Int, col: Int, remappingFunction: (T) -> T) {
        columns[col].mapAssign(row, remappingFunction)
    }

    override fun times(X: SparseVector<T>): SparseVector<T> {
        val result =  columns.first().new()
        for(j in 0 until nCols) {
            result.weightedPlusAssign(columns[j], X[j])
        }
        return result
    }

    // number of nonZeroEntries / (nRows*nCols)
    fun sparsity(): Double {
        val nonZeroEntriesSize = columns.sumBy { it.nonZeroEntries.size }
        return nonZeroEntriesSize.toDouble()/(nRows*nCols)
    }

    class CEntry<T>(override val col: Int, val entry: RemovableEntry<MutableMap.MutableEntry<Int,T>>, val isZero: T.()->Boolean): SparseMatrix.Entry<T> {
        override val row: Int
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

        override fun toString(): String {
            return "($row,$col)=$value"
        }
    }

//    class Entry<T>(override val col: Int, val colEntry: MutableMap.MutableEntry<Int,T>): SparseMatrix.Entry<T> {
//        override val row: Int
//            get() = colEntry.key
//        override val value: T
//            get() = colEntry.value
//        override fun setValue(newValue: T): T = colEntry.setValue(newValue) // TODO: Deal with set to zero
//    }


    // Iterator through a flattened list of MutableSparseVectors
//    class EntryIterator<T: Any>(
//        val vectors: List<MutableSparseVector<T>>,
//        var vectorIndex: Int,
//        var vectorIterator:MutableIterator<MutableMap.MutableEntry<Int,T>>): MutableIterator<Entry<T>> {
//
//        constructor(columns: List<MutableSparseVector<T>>): this(columns,0, columns[0].nonZeroEntries.entries.iterator())
//
//        override fun hasNext(): Boolean {
//            if(vectorIterator.hasNext()) return true
//            for(i in vectorIndex + 1 until vectors.size) {
//                if(!vectors[i].nonZeroEntries.isEmpty()) return true
//            }
//            return false
//        }
//
//        override fun next(): Entry<T> {
//            while(!vectorIterator.hasNext() && vectorIndex != vectors.lastIndex)
//                vectorIterator = vectors[++vectorIndex].nonZeroEntries.entries.iterator()
//            return Entry(vectorIndex, vectorIterator.next())
//        }
//
//        override fun remove() { vectorIterator.remove() }
//
//    }

}