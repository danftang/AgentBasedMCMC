package lib.sparseMatrix

import lib.vector.MutableSparseVector
import lib.vector.SparseVector

interface SparseRowMatrix<T: Any>: SparseMatrix<T> {
    val rows: List<MutableSparseVector<T>>

    override val nRows: Int
        get() = rows.size

    override val nonZeroEntries: MutableIterable<SparseMatrix.Entry<T>>
        get() = object: MutableIterable<Entry<T>> { override fun iterator() = EntryIterator(rows) }

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

//    fun clearRow(i: Int)
//    fun swapRows(i1: Int, i2: Int)
//
//    fun removeRows(rowsToRemove: Iterable<Int>)
//
//    fun replaceNonZeroElementsInRow(col: Int, map: (row: Int, value: Int) -> Int)
//
//    // row[i1] = row[i1] + weight*row[i2]
//    fun weightedRowPlusAssign(i1: Int, i2: Int, weight: Int) {
//        for(entry in rows[i2]) {
//            plusAssign(i1, entry.key, weight*entry.value)
//        }
//    }

    class Entry<T>(override val row: Int, val rowEntry: MutableMap.MutableEntry<Int,T>): SparseMatrix.Entry<T> {
        override val col: Int
            get() = rowEntry.key
        override val value: T
            get() = rowEntry.value
        override fun setValue(newValue: T): T = rowEntry.setValue(newValue)
    }

    class EntryIterator<T: Any>(
        val vectors: List<MutableSparseVector<T>>,
        var vectorIndex: Int,
        var vectorIterator:MutableIterator<MutableMap.MutableEntry<Int,T>>): MutableIterator<Entry<T>> {

        constructor(rows: List<MutableSparseVector<T>>): this(rows,0, rows[0].nonZeroEntries.entries.iterator()) {
            advanceToNextEntry()
        }

        override fun hasNext() = vectorIterator.hasNext()

        override fun next(): Entry<T> {
            val entry = vectorIterator.next()
            advanceToNextEntry()
            return Entry(vectorIndex, entry)
        }

        override fun remove() { vectorIterator.remove() }

        private inline fun advanceToNextEntry() {
            while(!vectorIterator.hasNext() && vectorIndex != vectors.lastIndex)
                vectorIterator = vectors[++vectorIndex].nonZeroEntries.entries.iterator()
        }
    }

}