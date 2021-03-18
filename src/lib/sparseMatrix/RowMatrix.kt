package lib.sparseMatrix

import lib.sparseVector.SparseVector

interface RowMatrix<T> {
    val rows: List<SparseVector<T>>
    val nCols: Int
    val nRows: Int
        get() = rows.size

}


operator fun <T> RowMatrix<T>.times(X: SparseVector<T>): SparseVector<T> {
    val result = X.new()
    for(i in rows.indices) {
        result[i] = rows[i].dotProduct(X)
    }
    return result
}