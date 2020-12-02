package lib.sparseMatrix

import lib.vector.MutableSparseVector
import lib.vector.SparseVector

interface SparseColMatrix<T: Any>: SparseMatrix<T> {
    val columns: List<MutableSparseVector<T>>

    override val nCols: Int
        get() = columns.size

    override val entries: Sequence<SparseMatrix.Entry<T>>
        get() = columns
            .asSequence()
            .mapIndexed { j, col ->
                col.nonZeroEntries.asSequence().map { SparseMatrix.Entry(it.key, j, it.value) }
            }
            .flatten()


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
        val result = columns.first().new()
        for(j in 0 until nCols) {
            result.weightedPlusAssign(columns[j], X[j])
        }
        return result
    }
}