package lib.sparseMatrix

import lib.vector.MutableSparseVector
import lib.vector.SparseVector

interface SparseColMatrix<T: Any>: SparseMatrix<T> {
    val columns: List<MutableSparseVector<T>>

    override val nCols: Int
        get() = columns.size

    override val nonZeroEntries: MutableIterable<SparseMatrix.Entry<T>>
        get() = object: MutableIterable<Entry<T>> { override fun iterator() = EntryIterator(columns) }

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

    class Entry<T>(override val col: Int, val colEntry: MutableMap.MutableEntry<Int,T>): SparseMatrix.Entry<T> {
        override val row: Int
            get() = colEntry.key
        override val value: T
            get() = colEntry.value
        override fun setValue(newValue: T): T = colEntry.setValue(newValue)
    }

    class EntryIterator<T: Any>(
        val vectors: List<MutableSparseVector<T>>,
        var vectorIndex: Int,
        var vectorIterator:MutableIterator<MutableMap.MutableEntry<Int,T>>): MutableIterator<Entry<T>> {

        constructor(columns: List<MutableSparseVector<T>>): this(columns,0, columns[0].nonZeroEntries.entries.iterator()) {
            advanceToNextEntry()
        }

        override fun hasNext() = vectorIterator.hasNext()

        override fun next(): Entry<T> {
            val sparseVecEntry = vectorIterator.next()
            val matrixEntry = Entry(vectorIndex, sparseVecEntry)
            advanceToNextEntry()
            return matrixEntry
        }

        override fun remove() { vectorIterator.remove() }

        private inline fun advanceToNextEntry() {
            while(!vectorIterator.hasNext() && vectorIndex != vectors.lastIndex)
                vectorIterator = vectors[++vectorIndex].nonZeroEntries.entries.iterator()
        }
    }

}