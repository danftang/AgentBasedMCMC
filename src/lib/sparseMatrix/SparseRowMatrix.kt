package lib.sparseMatrix

import lib.vector.MutableSparseVector
import lib.vector.SparseVector

interface SparseRowMatrix<T: Any>: SparseMatrix<T> {
    val rows: List<MutableSparseVector<T>>

    override val nRows: Int
        get() = rows.size

    override val entries: Sequence<SparseMatrix.Entry<T>>
        get() = rows
            .asSequence()
            .mapIndexed { i, row ->
                row.nonZeroEntries.asSequence().map { SparseMatrix.Entry(i, it.key, it.value) }
            }
            .flatten()

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

}